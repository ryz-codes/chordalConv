function [U] = minrank(X, clique, parent)
%MINRANK Minimum rank chordal completion

% Generate a postordering via recursive DFS
% TODO: replace with nonrecursive DFS!
[ch, root] = p2ch(parent); % Get children pointers
nr = numel(root);
post = cell(1,nr);
for i = 1:nr
    post{i} = dfs(root(i));
end
post = [post{:}];
function [ord] = dfs(v) % Recursive depth-first search
    if isempty(ch{v}) % Leaf
        ord = v;
    else
        % Traverse children
        ord = cell(1,numel(ch{v}));
        for vc = 1:numel(ch{v})
            ord{vc} = dfs(ch{v}(vc));
        end
        ord = [ord{:}, v];
    end
end

% Parse lengths
n = length(X);
omegas = cellfun(@numel, clique);
omega = max(omegas);

% Initialize
U = zeros(n,omega);
for i = post(end:-1:1) % reverse topological order
    J = clique{i}; par = parent(i);
    
    % Extract the submatrix for this clique
    Ui = zeros(omegas(i),omega);
    Ui(:,1:omegas(i)) = fact(X(J,J));
    if par == 0 % root, no overlaps
        U(clique{i},:) = Ui;
    else % Must figure out rotation
        % Parent clique
        Jp = clique{par};
        mask = ismember(J,Jp); maskp = ismember(Jp,J);
        
        % Overlapping rows
        Upover = U(Jp(maskp),:); Uover = Ui(mask,:);
        
        % Figure out rotation
        R = polar(Uover'*Upover);
        U(J(~mask),:) = Ui(~mask,:)*R;
    end
end

end
function U = fact(X)
[U,S,~] = eig(full(X+X'));
S = max(real(S),0);
U = U*sqrt(S)/sqrt(2);
end
function R = polar(Mat)
[U,~,V] = svd(Mat); R = U*V';
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
