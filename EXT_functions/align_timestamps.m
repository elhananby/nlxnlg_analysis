function [offset,amatch,bmatch] = align_timestamps(a,b,thresh)
%
% align_timestamps
%
% Align vectors of timestamps measured by two different clocks, using the 
% Needleman-Wunsch algorithm.
%
% Inputs:
%       a,b:    Vectors of timestamps, not necessarily the same length
%       thresh: Difference threshold
%
% Outputs:
%       offset: polyfit structure, the mean offset between amatch and
%           bmatch will be offset.mu(1);
%       amatch: matching timestamps from a.
%       bmatch: matching timestamps from b, same size as amatch.
%
% bmatch ~= amatch + polyval(offset.p,amatch,offset.S,offset.mu);
%
% a and b are timestamps, a subset of which denote the same events measured
% by clocks A and B, not synchronized with each other, and may drift linearly 
% with respect to each other. There can be spurious events or imprecise 
% timestamping. The algorithm compares between timestamp difference vectors. 
% If the difference vectors are within a specified threshold of each other, 
% it is considered a match. Constructs a scoring matrix and a direction 
% matrix for the NW algorithm, and then computes a maximum score path 
% through the scoring matrix. The path is decoded to form the optimal 
% alignment. Matching timestamps are fitted by a line to reveal the offset 
% between the clocks, and the clock drift. Timestamps from A can now be 
% shifted to approximate timestamps in B as:
%
% bTs = aTs + polyval(offset.p,aTs,offset.S,offset.mu);
%
% References: 
% http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
% http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html
% http://www.hrbc-genomics.net/training/bcd/Curric/PrwAli/node3.html
%
% Manu S. Madhav
% Written,  21-Apr-2015
% Modified, 02-Sep-2015, V1.0

% Vectorize
a = a(:)';
b = b(:)';

% Differences
da = diff(a);
db = diff(b);

m = length(da);
n = length(db);

gap = -1;
match = 1;
miss = -1;

% Construct scoring and direction matrices
S = zeros(m+1,n+1);
D = cell(m+1,n+1);
S(1,:) = [0,gap*(1:n)];
S(:,1) = [0;gap*(1:m)'];
[D{1,2:end}] = deal(3);
[D{2:end,1}] = deal(2);
for j = 2:m+1
    for k = 2:n+1
        if abs(db(k-1)-da(j-1))<=thresh
            val=match;
        else
            val=miss;
        end
        vec = [S(j-1,k-1) + val, S(j-1,k) + gap, S(j,k-1) + gap];
        S(j,k) = max(vec);
        D{j,k} = find(vec==S(j,k));
    end
end

% Find path through scoring grid
path = [m+1,n+1];
flag = 1;
while flag
    pathend = path(end,:);
    if pathend(1)==0 || pathend(2)==0
        continue;
    end
    d = D{pathend(1),pathend(2)};
    for k = 1:length(d)
        if d(1)==1
            newpathpoint = [pathend(1)-1,pathend(2)-1];
        elseif d(k)==2
            newpathpoint = [pathend(1)-1,pathend(2)];
        else
            newpathpoint = [pathend(1),pathend(2)-1];
        end

        path = [path;newpathpoint];   
    end

    if path(end,1)==1 && path(end,2)==1
        flag = 0;
    end
end
path = flipud(path);

% Matches only where there is a diagonal path and where the timestamps match
diffpath = diff(path);
matchidx = find(diffpath(:,1) & diffpath(:,2));

matchidx = matchidx(abs(da(path(matchidx+1,1)-1)-db(path(matchidx+1,2)-1))<=thresh);
matchidx = unique([matchidx;matchidx+1]);
matchidx = path(matchidx,:);

amatch = a(matchidx(:,1));
bmatch = b(matchidx(:,2));
diffvec = bmatch-amatch;

% Fit line to drift vector to form offset structure
[offset.p,offset.S,offset.mu] = polyfit(amatch,diffvec,1);
