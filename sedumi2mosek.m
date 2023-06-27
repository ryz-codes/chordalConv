function prob = sedumi2mosek(At, b, c, K)
% Convert SDP from SeDuMi format to MOSEK format
% Calling sequence 
%    [r,res] = mosekopt('minimize info',sedumi2mosek(A, b, c, K));
%    [x,y,rp,rd] = moseksol2sedumi(res, A, b, c, K); 
% is equivalent to
%    [x,y] = sedumi(A,b,c,K);
% Completely rewritten the code to make it run fast; added input checks
% Author: Richard Y Zhang <ryz@illinois.edu>
% Date:   May 20th, 2021

% This program is licenced under the BSD 2-Clause licence,
% contained in the LICENCE file in this directory.

% Input check: Cone K
if ~isfield(K,'f'), K.f = 0; end
if ~isfield(K,'l'), K.l = 0; end
if ~isfield(K,'q'), K.q = 0; end
if ~isfield(K,'s'), K.s = 0; end

% Input check: Matrix At
order = K.f+K.l+sum(K.q)+sum(K.s.^2);
if size(At,2) == order
    % A is provided instead of At, so take transpose.
    At = At';
elseif size(At,1) ~= order
    % Dimensions do not match, the code is going to error
    error('The provided At matrix does not have the right number of columns compared to the order of the cone K');
end

% Input check: Vector b
b = b(:);
if numel(b) ~= size(At,2)
    error('The provided b vector does not have the right number of elements compared to the number of columns in the matrix At');
end

% Input check: Vector c
c = c(:);
if numel(c) ~= order
    error('The provided c vector does not have the right number of elements compared to the order of the cone K');
end

% Number of rows in At occupied by the f,l,q blocks
flq = K.f+K.l+sum(K.q);

% Form f, l, q blocks in a simple way
prob.a = At(1:flq,:)';
prob.c = c(1:flq);
prob.blx = [-inf(K.f,1);zeros(K.l,1);-inf(sum(K.q),1)]; % Primal lower-bound
prob.bux = inf(flq,1); % Primal upper-bound
prob.blc = b; % Linear constraint lower-bound
prob.buc = b; % Linear constraint upper-bound
if K.q(1)>0    
    prob.cones.type = zeros(1,length(K.q));
    prob.cones.subptr = cumsum([1, K.q(1:end-1)]);
    prob.cones.sub = K.f + K.l + (1:sum(K.q));
end

% Form SDP blocks 
prob.bardim = K.s;
prob.bara = splitSDPblocks(At(flq+1:end,:), K.s);
prob.barc = splitSDPblocks(c(flq+1:end,:), K.s);
end

function bar = splitSDPblocks(At, s)
% Put sedumi style At matrix SDP blocks into MOSEK format
% inputs:   At     - sedumi At matrix
%           s      - the orders of the SDP cones (this is K.s)
% MOSEK format requires the following value pairs
%       A_{i,j} [k,l] = val
% in order to implement the following constraints
%       sum_j <A_{i,j}, X_j> for all i
% where each X_j is in cone(j)
[allcol, allrow, allcon, allvar, allval] ...
    = deal(cell(length(s), 1)); % Preallocate

% Convert into i,j,k format, and then 
[ii,jj,kk] = find(At);
[N,m] = size(At);

% reorder from column-raster into row-raster
[ii,perm] = sort(ii);
jj = jj(perm); kk = kk(perm);

% Isolate blocks according to what row they are in
rowpnt = 1; % Pointer to the current row being looked at
valpnt = 1; % Pointer to the current nonzero being looked at
for k = 1:length(s)  
    % Decompose the k-th SDP block of A
    %   A{k}*x{k} = [<A{k,j}, X{k}>] where X{k} in SDP(K.s(k))
    blocklen = s(k)^2;
    
    % Efficiently implement the following logic
    %  [this_ii2,this_jj2,this_kk2] = find(At(rowpnt:rowpnt + blocklen - 1,:));
    
    % Want to find all nonzeros with rowpnt <= i(t) < newrowpnt
    newrowpnt = rowpnt + blocklen; 
    found = false; 
    for newvalpnt = valpnt:length(ii)
        if ii(newvalpnt) >= newrowpnt
            found=true; break; 
        end
    end
    if ~found, newvalpnt = length(ii)+1; end
    
    % The nonzero at t = newvalpnt satisfies i(t) = newrowpnt
    this_ii = ii(valpnt:newvalpnt-1) - rowpnt + 1;
    this_jj = jj(valpnt:newvalpnt-1);
    this_kk = kk(valpnt:newvalpnt-1);
    
    % Increment pointers
    valpnt = newvalpnt;
    rowpnt = newrowpnt;

    % Map indices in R^{n^2} back to row/column subscripts in R^{n x n}
    [row,col] = ind2sub([s(k), s(k)], this_ii); 
    mask = row >= col; % lower-triangular mask
    
    % Store into preallocated storage
    allcol{k} = col(mask); % column subscripts
    allrow{k} = row(mask); % row subscripts
    allval{k} = this_kk(mask);  % values
    allvar{k} = this_jj(mask);  % which constraint / dual variable?
    allcon{k} = k*ones(sum(mask),1); % which SDP block?
end
bar.subi = cat(1,allvar{:});
bar.subj = cat(1,allcon{:});
bar.subk = cat(1,allrow{:});
bar.subl = cat(1,allcol{:});
bar.val = cat(1,allval{:});
end

function [this_ii,this_jj,this_kk] = findblk(blkpnt,blocklen,ii,jj,kk,N,m)
% Try to optimize the following line
%   [this_ii,this_jj,this_kk] = find(At(blkpnt:blkpnt + blocklen - 1,:));
% where At = sparse(ii,jj,kk,n2,m)
    ii = ii - blkpnt;

    % The following code is an optimized version of the following
      blkstart = find(ii >= 1, 1, 'first'); 
      blkend = find(ii >= blocklen+1, 1, 'first'); 
      if isempty(blkend) % Find returns empty if it reaches the end
          blkend = length(ii)+1;
	  end
%     found = false; thresh = blkpnt + blocklen;
%     for nxtblkpnt = blkpnt:length(ii)
%         if ii(nxtblkpnt) >= thresh
%             found=true; break; 
%         end
%     end
%     if ~found, nxtblkpnt = length(ii)+1; end
    
    % Isolate the k-th block-row of At.
    this_ii = ii(blkstart:blkend-1);
    this_jj = jj(blkstart:blkend-1);
    this_kk = kk(blkstart:blkend-1);
end

% function bar = splitSDPblocks(At, s, offset)
% % Put sedumi style At matrix SDP blocks into MOSEK format
% % inputs:   At     - sedumi At matrix
% %           offset - number of rows to skip in At
% %           s      - the orders of the SDP cones (this is K.s)
% if nargin < 3 || isempty(offset)
%     pointer = 0;
% else
%     pointer = offset;
% end
% [allcol, allrow, allcon, allvar, allval] ...
%     = deal(cell(length(s), 1)); % Preallocate
% for k = 1:length(s)  
%     % Decompose the k-th SDP block of A
%     %   A{k}*x{k} = [<A{k,j}, X{k}>] where X{k} in SDP(K.s(k))
%     blocklen = s(k)^2;
%     [ii,jj,kk] = find(At(pointer+(1:blocklen),:));
%     pointer = pointer + blocklen;
%     
%     % Map indices in R^{n^2} back to row/column subscripts in R^{n x n}
%     [row,col] = ind2sub([s(k), s(k)], ii); 
%     mask = row >= col; % lower-triangular mask
%     
%     % Store into preallocated storage
%     allcol{k} = col(mask); % column subscripts
%     allrow{k} = row(mask); % row subscripts
%     allval{k} = kk(mask);  % values
%     allvar{k} = jj(mask);  % which constraint / dual variable?
%     allcon{k} = k*ones(sum(mask),1); % which SDP block?
% end
% bar.subi = cat(1,allvar{:});
% bar.subj = cat(1,allcon{:});
% bar.subk = cat(1,allrow{:});
% bar.subl = cat(1,allcol{:});
% bar.val = cat(1,allval{:});
% end
