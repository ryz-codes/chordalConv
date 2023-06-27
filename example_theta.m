n = 1e4;  % Number of vertices in the graph
k = 34;   % Treewidth of the graph
sp = 1.5; % Ratio of edges to vertices

% Generate a random graph with specified treewidth
Adj = kTree(n,k,sp/k); 

% Generate an instance of the Lovasz Theta problem on this graph
% The data A,b,c,K are given in SeDuMi format
[A,b,c,K] = genTheta(Adj); lb = b; ub = b;

% Solve the SDP
[U,y,info] = solveChordConv(c,A,lb,ub);