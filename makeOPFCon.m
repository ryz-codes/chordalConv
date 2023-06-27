function [At, lb, ub, c, xsol] = makeOPFCon(data,opts)
%MAKEOPFCON Given MATPOWER case file, make A, b, C matrices so that the
%standard state estimation problem can be posed as the semidefinite program
%   min   c' * x
%   s.t.  lb <= A' * x  <= ub
%         mat(x) >= 0 
% MATPOWER must be in the search path.

% OPF formulation:
%   Loads are fixed PQ
%   All bus voltages constrained from 0.95 to 1.05
%   Generators subject to the D curve. 
% minimize generation cost, proportional to P (and oblivious to Q)

% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   May 20th, 2023

% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in this directory.


% Default values
VOLT_TOL = 0.05;
LOAD_TOL = 0.2;
MARGIN   = 0.2;
LINE_TOL = 1;

% Load options
if nargin < 2 || isempty(opts)
    opts = struct;
end
doreal = isfield(opts,'doreal') && opts.doreal;
if isfield(opts,'volt_tol') && ~isempty(opts.volt_tol)
    VOLT_TOL = opts.volt_tol; end
if isfield(opts,'load_tol') && ~isempty(opts.load_tol)
    LOAD_TOL = opts.load_tol; end
if isfield(opts,'margin') && ~isempty(opts.margin)
    MARGIN = opts.margin; end
if isfield(opts,'line_tol') && ~isempty(opts.line_tol)
    LINE_TOL = opts.line_tol; end


% Input
%    data  -- MATPOWER case file
% Output
%  A,b,c,K -- Semidefinite program in SeDuMi format
%             The first 2*nbus constraints in A and b correspond to the
%             usual power flow problem. The remaining constraints are the
%             remaining state estimation measurements
%    xsol  -- Unique rank-1 solution
%             for the powerflow problem

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

% Presolve powerflow.
data = runpf(data);

%==========================================================================
% Standard matpower commands
mpc = loadcase(data);
mpc = ext2int(mpc);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
%==========================================================================
    
% BUSES
[Ybus, Yf] = makeYbus(baseMVA, bus, branch);
Sbus = V0 .* conj(Ybus*V0); n = numel(V0); 

% Define the following matrix-valued operators at the j-th bus
%   E = @(j) Id(:,j)*Id(:,j)'
%   S = @(j) conj(Ybus)*E(j); 
%   P = @(j) (S(j) + S(j)')/2;
%   Q = @(j) (S(j) - S(j)')/2i;
% Then our code below implements the following in a vectorized way
%   Vbus(:,j) = vec(E(j));
%   Pbus(:,j) = vec(P(j));
%   Qbus(:,j) = vec(Q(j));
[ii,jj,kk] = find(conj(Ybus));
if doreal
    Vbus = sparse(1:2*n+1:(4*n^2), [1:n,1:n], 1, 4*n^2, n);  
    Pbus = sparse([ii + (jj-1)*2*n,     (ii-1)*2*n + jj, ...      (1,1) block
                   ii+n + (jj+n-1)*2*n, (ii+n-1)*2*n + jj+n, ...  (2,2) block
                   ii+n + (jj-1)*2*n,   (ii+n-1)*2*n + jj, ...   (2,1) block
                   ii + (jj+n-1)*2*n,   (ii-1)*2*n + jj+n], ...  (1,2) block
                   jj*ones(1,8), ...
                   [real(kk)*[1,1,1,1], imag(kk)*[1,1,-1,-1]], 4*n^2, n);
    Qbus = sparse([ii + (jj-1)*2*n,     (ii-1)*2*n + jj, ...      (1,1) block
                   ii+n + (jj+n-1)*2*n, (ii+n-1)*2*n + jj+n, ...  (2,2) block
                   ii+n + (jj-1)*2*n,   (ii+n-1)*2*n + jj, ...   (2,1) block
                   ii + (jj+n-1)*2*n,   (ii-1)*2*n + jj+n], ...  (1,2) block
                   jj*ones(1,8), ...
                   [imag(kk)*[1,1,1,1],real(kk)*[-1,-1,1,1]], 4*n^2, n);
else
    Vbus = sparse(1:n+1:n^2,1:n,1,n^2,n);  % = vec(ej*ej') at j-th col
    Pbus = sparse([ii + (jj-1)*n, (ii-1)*n + jj], jj*[1,1], [kk,conj(kk)]/2, n^2, n);
    Qbus = sparse([ii + (jj-1)*n, (ii-1)*n + jj], jj*[1,1], [kk,-conj(kk)]/2i, n^2, n);
end

% Voltage constraints
Av = Vbus;
lbv = abs(V0(:)).^2*(1-VOLT_TOL)^2;
ubv = abs(V0(:)).^2*(1+VOLT_TOL)^2;

% PQ constraints
pq = bus(:,BUS_TYPE) == PQ;
Apq = [Pbus(:,pq),  Qbus(:,pq)]; 
bpq = [real(Sbus(pq)); imag(Sbus(pq))];
ubpq = bpq + abs(bpq)*LOAD_TOL;
lbpq = bpq - abs(bpq)*LOAD_TOL;

% Generators
pv = bus(:,BUS_TYPE) == PV;
Apv  = [Pbus(:,pv), Qbus(:,pv)]; 
Pmax = abs(real(Sbus(pv)))*(1+MARGIN);
Qmax = abs(imag(Sbus(pv)))*(1+MARGIN);
ubpv = [Pmax;      Qmax];
lbpv = [Pmax*0;   -Qmax];

% Cost function
genid = gen(:,GEN_BUS);
c = sum(Pbus(:,genid),2);

%==========================================================================
% BRANCHES
% Define the following matrix-valued operators at the j-th bus
%   E = @(j) Id(:,j)*Id(:,j)'
%   S = @(j) Id(:,from(j)) * conj(Yf(j,:)); 
%   P = @(j) (S(j) + S(j)')/2;
%   Q = @(j) (S(j) - S(j)')/2i;
% Then our code below implements the following in a vectorized way
%   Pbranch(:,j) = vec(P(j));
%   Qbus(:,j)  = vec(Q(j));
[ii,jj,kk] = find(Yf');
from = branch(:,F_BUS); m = numel(from);
fjj = from(jj);
if doreal
    Pbranch = sparse([ii + (fjj-1)*2*n,     (ii-1)*2*n + fjj, ...      (1,1) block
                   ii+n + (fjj+n-1)*2*n, (ii+n-1)*2*n + fjj+n, ...  (2,2) block
                   ii+n + (fjj-1)*2*n,   (ii+n-1)*2*n + fjj, ...   (2,1) block
                   ii + (fjj+n-1)*2*n,   (ii-1)*2*n + fjj+n], ...  (1,2) block
                   jj*ones(1,8), ...
                   [real(kk)*[1,1,1,1], imag(kk)*[1,1,-1,-1]], 4*n^2, m);
else
    Pbranch = sparse([ii + (fjj-1)*n, (ii-1)*n + fjj], jj*[1,1], [kk,conj(kk)]/2, n^2, m);
    %Qbranch = sparse([ii + (from(jj)-1)*n, (ii-1)*n + from(jj)], jj*[1,1], [kk,-conj(kk)]/2i, n^2, m);
end
Sline = V0(from) .* conj(Yf*V0);

% Line constraints (from)
if isfinite(LINE_TOL)
    Al = Pbranch; bl = real(Sline);
    ubl = +abs(bl)*(1+LINE_TOL);
    lbl = -abs(bl)*(1+LINE_TOL);
else
    Al = []; ubl = []; lbl = [];
end

%==========================================================================
% ASSEMBLE
At = [Av,  Apq, Apv, Al];
ub = [ubv; ubpq; ubpv; ubl];
lb = [lbv; lbpq; lbpv; lbl];
xsol = V0;
end




