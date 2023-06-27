function [Adj] = kTree(n,k,ratio)
% KTREE    Generates random partial k-tree.
%
% Returns adjacency matrix and cell array with cliques.

assert(k > 0 && round(k) == k,'k must be a positive integer')
assert(n > 0 && round(n) == n && n > k,'n must be a positive integer (n>k)')
if nargin < 3, ratio = 1;end

% Generate cliques
C = cell(n-k,1);
C{1} = 1:k+1;
for j = 2:n-k
    mask = true(1,k+1);
    mask(randi(k+1)) = false;
    C{j} = [C{randi(j-1)}(mask) k+j];
end

% Generate complete adjacency matrix
A = sparse(cat(1,C{:}),repmat(1:n-k,1,k+1),true,n,n-k);
Adj = logical(A*A');

% Subsample edges 
if ratio < 1
    [I,J] = find(tril(Adj,-1)); 
    m = numel(I);
    sel = rand(m,1) <= ratio;
    Adj = sparse(I(sel),J(sel),true,n,n);
    Adj = Adj | Adj';
end

% Apply permutation
p = n:-1:1;
Adj = Adj(p,p);
% for j = 1:n-k
%     C{j} = p(C{j});
% end

end

