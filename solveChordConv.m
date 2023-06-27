function [U, ys, info] = solveChordConv( c, A, lb, ub, opts)
%SOLVECHORDCONV Solve chordal conversion
% Reformulate the primal standard-form problem
%   min c'x s.t. lb <= A*x <= ub, x \in SDP(n)
% to the dual standard-form problem
%   max -c'y s.t. lb <= A*x <= ub, x[Jj,Jj] \in SDP(nj) for all j
% and solve the problem using MOSEK.

% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   May 20th, 2023

% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in this directory.

if nargin < 5 || isempty(opts)
    opts = struct;
end
if isfield(opts,'perm')
    perm = opts.perm;
else
    perm = [];
end
if isfield(opts,'verbose')
    verbose = opts.verbose;
else
    verbose = 1;
end
if isfield(opts,'norecover')
    norecover = true;
else
    norecover = false;
end

% Assemble data
tic; [Amat,bvec,cvec,Kcone,info] = chordConv( c, A, lb, ub, perm);
info.time.conv = toc;

% SeDuMi format to MOSEK format
tic; prob = sedumi2mosek(Amat,bvec,cvec,Kcone);
info.time.sed2mos = toc;

% Call MOSEK
param = struct;
if isfield(opts,'natural') && opts.natural
    param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
    param.MSK_IPAR_INTPNT_ORDER_METHOD = 'MSK_ORDER_METHOD_NONE';
end
if isfield(opts,'pfeas') && ~isempty(opts.pfeas)
    param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = pfeas;
end
param.MSK_IPAR_INTPNT_SOLVE_FORM = 'MSK_SOLVE_PRIMAL';
msk_cmd = 'minimize info';
if verbose ==0, msk_cmd = [msk_cmd, ' echo(0)']; end
tic;[r,res] = mosekopt(msk_cmd,prob,param);
info.time.mosek = toc;
info.msk.time_presolve = res.info.MSK_DINF_PRESOLVE_TIME;
info.msk.time_order = res.info.MSK_DINF_INTPNT_ORDER_TIME;
info.msk.time_ipm = res.info.MSK_DINF_INTPNT_TIME;
info.msk.iter_ipm = res.info.MSK_IINF_INTPNT_ITER;
info.msk.time_total = res.info.MSK_DINF_OPTIMIZER_TIME;
info.msk.info = res.info;

% MOSEK format to SeDuMI format
tic; [x,y,rp,rd,cx,by] = moseksol2sedumi(res, Amat,bvec,cvec,Kcone); 
info.time.mos2sed = toc;
pinf = rp / (1+norm(bvec)); % |Ax - b| / (1 + |b|)
dinf = rd / (1+norm(cvec)); % [Ay - c]_+ / (1 + |c|)
dgap = abs(cx-by)/(1+abs(cx)+abs(by));
digits  = -log10(pinf + dinf + dgap);

if norecover % Skip recovering U and y
    U = []; ys = [];
    optX = by; optU = by; opty = cx;
    info.time.recov_x = 0;
    info.time.recov_u = 0;
    info.time.recov_y = 0;
else   
    % Recover sparse primal variable
    tic; [Ei,Ej] = find(info.L); 
    n = info.n; nr = numel(Ei); 
    if info.doreal
        Xs = sparse(Ei, Ej, y, n, n);    
    else
        y(nr+1:end) = 1j*y(nr+1:end); s = Ei>Ej; 
        Xs = sparse([Ei;Ei(s)], [Ej;Ej(s)], y, n, n);
    end
    Xs = Xs + Xs';
    info.time.recov_x = toc;

    % Recover chordal low-rank completion
    tic; U = minrank(Xs,info.clique,info.parent);
    info.time.recov_u = toc;
    U(info.perm,:) = U;
    Xs(info.perm,info.perm) = Xs;

    % Recover dual variables
    tic; m = size(A,2); ys = zeros(m,2);
    do_lb = isfinite(lb); do_ub = isfinite(ub);
    ys(:,1) = -x(1:nnz(do_lb));
    ys(:,2) = -x(nnz(do_lb)+(1:nnz(do_ub)));
    ys = info.D*ys;
    info.time.recov_y = toc;

    % The three optimal values
    optX = full(c(:)'*Xs(:)); n = size(U,1);
    optU = U(:)'*reshape(reshape(c,n,n)*U,[],1);
    opty = ys(:)'*[-lb;ub];
end

info.sol.optX = optX;
info.sol.optU = optU;
info.sol.opty = opty;
info.sol.pinf = pinf;
info.sol.dinf = dinf;
info.sol.dgap = dgap;
info.sol.digits = digits;

if verbose > 0
fprintf('prim(X):      %1.12e\n',optX);
fprintf('prim(U*U''):   %1.12e\n',optU);
fprintf('dual(y):      %1.12e\n',opty);
fprintf('-----------------------\n');
fprintf('lb<=A(U*U'')<=ub:     %1.12e\n',pinf);
fprintf('c - A''(y)>=0:        %1.12e\n',dinf);
fprintf('prim(U*U'')=dual(y):  %1.12e\n',dgap);
fprintf('Accurate digits:     %1.12e\n',digits);
fprintf('-----------------------\n');
fprintf('chordConv:    %1.6e sec\n',info.time.conv);
fprintf('sedumi2mosek: %1.6e sec\n',info.time.sed2mos);
fprintf('mosek:        %1.6e sec\n',info.time.mosek);
fprintf('   presolve   %1.6e sec\n',info.msk.time_presolve);
fprintf('   order      %1.6e sec\n',info.msk.time_order);
fprintf('   ipm        %1.6e sec\n',info.msk.time_ipm);
fprintf('   total      %1.6e sec\n',info.msk.time_total);
fprintf('mosek2sedumi: %1.6e sec\n',info.time.mos2sed);
fprintf('recover X:    %1.6e sec\n',info.time.recov_x);
fprintf('recover U:    %1.6e sec\n',info.time.recov_u);
fprintf('recover y:    %1.6e sec\n',info.time.recov_y);
fprintf('-----------------------\n');
perit = info.msk.time_ipm/info.msk.iter_ipm;
fprintf('Iterations:   %d \n',info.msk.iter_ipm);
fprintf('Per iter:     %1.6e sec\n',perit);
end
end
