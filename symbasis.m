function V = symbasis(pattern, perm, doreal)
%SYMBASIS Basis for the set of real symmetric matrices within the linear 
% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   July 28th, 2018

% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in this directory.

if isscalar(pattern) % input a scalar for dense basis
    n = pattern;
    pattern = sparse(tril(ones(n)));
else % input a sparsity pattern
    n = length(pattern);
    pattern = tril(pattern+pattern');
    pattern(1:n+1:end) = 1;
end
if nargin < 2 || isempty(perm)
    perm = 1:n;
end
if nargin < 3 || isempty(doreal)
    doreal = true;
end

% Enumerate the subscripts of a lower triangular matrix
% space of all n x n matrices.
ndof = nnz(pattern);
[Ei,Ej] = find(pattern);

% Convert the subscrips into indices
lower = perm(Ei) + n*(perm(Ej)-1);
upper = n*(perm(Ei)-1) + perm(Ej);

% Number according to DOFs
jj = (1:ndof)';

% Scaling to make the matrix unitary
kk = sqrt(0.5)*ones(ndof,2);
kk(lower==upper,:) = 0.5*ones(n,2);

% Output
V = sparse([lower,upper],[jj,jj],kk,n^2,ndof);

% Complex basis
if nargin > 2 && ~doreal
    pattern_im = tril(pattern,-1);
    ndof = nnz(pattern_im); %Skew symmetric part

    % Enumerate the subscripts of a lower triangular matrix
    [Ei,Ej] = find(pattern_im);

    % Convert the subscrips into indices
    lower = perm(Ei) + n*(perm(Ej)-1);
    upper = n*(perm(Ei)-1) + perm(Ej);

    % Number according to DOFs
    jj = (1:ndof)';

    % Scaling to make the matrix unitary
    kk = 1i*[sqrt(0.5)*ones(ndof,1);
             -sqrt(0.5)*ones(ndof,1)];
    
    % Antihermitian part
    V = [V,sparse([lower,upper],[jj,jj],kk,n^2,ndof)];
end
end