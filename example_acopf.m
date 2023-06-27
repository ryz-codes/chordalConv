% NOTE: MATPOWER (matpower.org) should already be installed.

%mpc = case3120sp;
mpc = case13659pegase;
%mpc = case_ACTIVSg25k;
%mpc = case_ACTIVSg70k;

% Generate acopf problem
opts = struct('line_tol',inf,'load_tol',0.05);
[A,lb,ub,c] = makeOPFCon(mpc,opts);

% Solve the SDP
[U,y,info] = solveChordConv(c,A,lb,ub);