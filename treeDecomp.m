function [L, clique, parent] = treeDecomp( Adj)
%TREEDECOMP Given graph adjacency matrix, compute tree decomposition by
%eliminating with canonical ordering

% Input:
%   Adj     - Graph adjacency matrix

% Output:
% structure T with fields:
%   L       - Filled lower-triangular structure
%   cliques - maximal cliques, ordered in topological ordering
%   parent  - parent structure for the supernodal elimination tree

% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   Feb 27, 2018
% Reference: 
% [VA] Vandenberghe, L., & Andersen, M. S. (2015). Chordal graphs and 
% semidefinite optimization. Foundations and Trends in Optimization, 
% 1(4), 241-433.
% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in the home directory.

% Input checks
n = size(Adj,1);
assert(size(Adj,2) == n, 'Adjacency matrix must be square');
Adj = logical(Adj); % convert to logical
Adj = Adj | Adj'; % symmetricize
Adj(1:n+1:end) = true; % fill the diagonal

% Compute elimination tree [VA, p.41] which is also a clique tree
% decomposition. 
[~,~,parent,~,L] = symbfact(Adj,'sym','lower');
clique = cell(1,n);
for i = 1:n
    clique{i} = find(L(:,i));
end

% Merge redundant cliques into supernodes.
%[clique, parent] = supernode(clique, parent);
end

function [clique2, parent2] = supernode(clique, parent)
%SUPERNODE Maximal supernodes and supernodal elimination tree
% Vandenberghe and Andersen Algorithm 4.1

% Input checks
n = numel(parent);
assert(numel(clique) == n, 'mismatch in number of elements');
assert(all(parent > (1:n) | parent ==0), 'cliques must be given in postordering');

% Data structures
deg = cellfun(@numel,clique);
ch = p2ch(parent);

isuper = zeros(1,n);
parent2 = zeros(1,n);
repre = zeros(1,n);
ell = 0;
for v = 1:n
    % Check to see if we should make a new supernode
    makeNew = true;
    for w = ch{v}
        if deg(w) == deg(v) + 1
            makeNew = false;
            break;
        end
    end
    % If yes, create a new supernode u. 
    % If not, get supernode u to add to
    if makeNew
        ell = ell+1;
        u = ell;
        repre(u) = v;
    else
        u = isuper(w);
    end
    % Add v to supernode u
    isuper(v) = u;
    for w = ch{v}
        z = isuper(w);
        if z ~= u
            parent2(z) = u;
        end
    end
end
parent2 = parent2(1:ell);
clique2 = clique(repre(1:ell));
end
function [ch, root] = p2ch(p)
% Convert the list of parent vectors into a cell array of children vectors.
n = length(p);
ch = cell(size(p));
root = [];
for i = 1:n
    par = p(i);
    if par > 0
        % Push child into parent
        ch{par} = [ch{par}, i];
    elseif par == 0
        root = [root, i];
    end
end
end

