function [Amat,bvec,cvec,Kcone,info] = chordConv( c, At, lb, ub, perm)
%CHORDCONV Chordal conversion
% Reformulate the primal standard-form problem
%   min c'x s.t. lb <= A*x <= ub, x \in SDP(n)
% to the dual standard-form problem
%   max -c'y s.t. lb <= A*x <= ub, x[Jj,Jj] \in SDP(nj) for all j
% and output in SeDuMi format.

% Author: Richard Y Zhang <ryz@illinois.edu>

% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in this directory.

verbose = 2;

% Defaults
if size(At,2) > size(At,1), At = At.'; end
if nargin < 4 || isempty(ub)
    ub = lb;
elseif isempty(lb)
    lb = ub;
end
if nargin < 5, perm = []; end
lb = lb(:); ub = ub(:); c = c(:);
assert(all(lb <= ub), 'All lower-bounds must be smaller than the upper-bounds!');

%--------------------------------------------------------------------------
% Tree decomposition
%--------------------------------------------------------------------------
if verbose > 0
	fprintf('Beginning tree decomposition....');
end

% Get adjacency pattern
Adj = any([At,c],2);
n = floor(sqrt(length(Adj)));
assert(length(Adj) == n^2, 'number of rows in A must be a perfect square');
Adj = reshape(Adj, n, n);

% Compute the tree decomposition
if isempty(perm)
    perm = amd(Adj); 
else
    assert(length(perm) == n,'perm is wrong length!!');
end
[L, clique, parent] = treeDecomp(Adj(perm,perm));
Ncliques = length(clique);
Ks = cellfun(@length,clique);
frontsize = max(Ks);
if verbose > 0
	fprintf('complete\n');
    fprintf(' Num. cliques: %d \n Max cliques size: %d\n', ...
        Ncliques, frontsize);
end

% Store tree decomposition info
info = struct;
info.Adj = Adj;
info.n = n;
info.perm = perm;
info.L = L; 
info.clique = clique; 
info.parent = parent;
info.frontsize = frontsize;

%--------------------------------------------------------------------------
% Form matrices
%--------------------------------------------------------------------------
% Basis to go between svecE and vec(R^{n x n})
[Ei,Ej] = find(L);
At = [At,c]; % Merge A and c for this step, and split after
if isreal(At)
    % Split conic Y[Jj,Jj] >=0 and then
    % project vec(S^n) down to svec_E(S)
    P = cell(1,Ncliques);
    for j = 1:Ncliques
        Pj = sparse(clique{j},1:Ks(j),1,n,Ks(j));
        P{j} = kron(Pj,Pj); 
    end
	% [Vsym] = symbasis(L);
	% Psym = Vsym'*Psym;
    Psym = [P{:}];
	lower = Ei + n*(Ej-1);
	upper = n*(Ei-1) + Ej;
	Psym = Psym(lower,:) + Psym(upper,:);

    % Project vec(S^n) down to svec_E(S)
	% [Vsym2] = symbasis(L, perm);
    % Asym = Vsym2'*A;
    % csym = Vsym2'*c;
	lower = perm(Ei) + n*(perm(Ej)-1);
	upper = n*(perm(Ei)-1) + perm(Ej);
	Asym = At(lower,:) + At(upper,:);
    info.doreal = true;
else    
    % Split conic Y[Jj,Jj] >=0 and then
    % project vec(S^n) down to svec_E(S)
    P = cell(1,Ncliques); R = cell(1,Ncliques);
    for j = 1:Ncliques
        Pj = sparse(clique{j},1:Ks(j),1,n,Ks(j));
        P{j} = kron(Pj,Pj);
        Vj = symbasis(Ks(j),[],false);
        Vj2 = realembed(Vj);
        R{j} = sparse(Vj*Vj2');
    end
	% [Vsym] = symbasis(L);
	% Psym = Vsym'*Psym;
	lower = Ei + n*(Ej-1);
	upper = n*(Ei-1) + Ej;
	s = Ei > Ej; % The strict off-diagonal
    Psym = [P{:}]; R = blkdiag(R{:});
	Plow = Psym(lower,:)*R; Pupp = Psym(upper,:)*R;
	Psym = [real(Plow + Pupp); 
		    imag(Plow(s,:) - Pupp(s,:))];
    
    
    % Project vec(S^n) down to svec_E(S)
    lower = perm(Ei) + n*(perm(Ej)-1);
	upper = n*(perm(Ei)-1) + perm(Ej);
	Alow = At(lower,:); Aupp = At(upper,:);
	Asym = [real(Alow + Aupp);
			imag(Alow(s,:) - Aupp(s,:))];
        
    Ks = 2*Ks; info.doreal = false;
end
% Split back up
csym = Asym(:,end); Asym = Asym(:,1:end-1);

%--------------------------------------------------------------------------
% Assemble and output
%--------------------------------------------------------------------------
% Apply diagonal preconditioning to AtA
D = diag(sparse(1./sqrt(sum(Asym.^2,1))));
info.D = D;

% Process infinite bounds
bvec = -csym(:);
do_lb = isfinite(lb); do_ub = isfinite(ub);
Amat = [-D*Asym(:,do_lb)'; D*Asym(:,do_ub)'; -Psym']; 
cvec = [-D*lb(do_lb); D*ub(do_ub); zeros(size(Psym,2),1)];
Kcone = struct;
Kcone.l = nnz(do_lb) + nnz(do_ub);
Kcone.q = 0;
Kcone.s = Ks;
   
end

function [V2] = realembed(V)
%Embedding the size-n hermitian PSD cone into the
%size-2n real symmetric PSD cone.
V = full(V);
[n2, m] = size(V);
V2 = zeros(4*n2, m);
n = floor(sqrt(n2));
for i = 1:m
    tmp = reshape(V(:,i),n,n);
    tmp2 = [real(tmp), -imag(tmp); imag(tmp), real(tmp)];
    V2(:,i) = tmp2(:);
end
V2 = sparse(V2);
end

