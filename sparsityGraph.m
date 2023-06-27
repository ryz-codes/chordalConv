function [Adj,Adjext] = sparsityGraph(c,At)
%SPARSITYGRAPH Adjacency matrix for the sparsity graph and extended 
% sparsity graph 
% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   May 20th, 2021
n = floor(sqrt(size(At,1)));
m = size(At,2);
assert(size(At,1) == n^2, 'number of rows in A must be a perfect square');

% Get adjacency pattern
Adj = reshape(any([c,At],2), n, n);
Adj = Adj | Adj';
Adj(1:n+1:end) = true;

if nargout>1

    % Construct the sparse vectorization basis for merging columns
    [Ei, Ej] = find(Adj); d = numel(Ei);
    blkstart = cumsum([1,sum(Adj,1)]); 
    [ii,jj] = deal(cell(1,n));
    for j = 1:n
        ii{j} = blkstart(j):blkstart(j+1)-1;
        jj{j} = j*ones(1,blkstart(j+1)-blkstart(j));
    end
    merge = sparse([ii{:}],[jj{:}],true,d,n,d);
    
    % Compute support of each constraint
    Apattern = At(sub2ind([n,n],Ei,Ej),:) ~=0;
    Supp = merge'*Apattern > 0;
    
    % Get extended adjacency pattern
    Adjext = Adj + Supp*Supp' > 0;
end
end

