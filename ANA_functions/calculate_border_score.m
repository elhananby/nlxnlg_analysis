function [highWallIdx, perc] = calculate_border_score(c)

spkN = length(c.timestamps);

% get distance for each wall
spkDist = [c.posx_c, 96.5 - c.posx_c, c.posy_c, 96.5 - c.posy_c];
spkDist(spkDist < 0) = NaN;

% find minimal distance for each spike and wall (ind)
% 1 - east; 2 - west; 3 - south; 4 - north
[minDist, minInd] = min(spkDist, [], 2);

% initialize distnace matrix
distVec = 1:1:96.5;

westWall = minDist(minInd == 1);
eastWall = minDist(minInd == 2);
southWall = minDist(minInd == 3);
northWall = minDist(minInd == 4);

perc = zeros(4, length(distVec));

for ii = distVec
    westDist = westWall(westWall <= ii);
    eastDist = eastWall(eastWall <= ii);
    southDist = southWall(southWall <= ii);
    northDist = northWall(northWall <= ii);
    
    perc(1, ii) = length(westDist)/spkN;
    perc(2, ii) = length(eastDist)/spkN;
    perc(3, ii) = length(southDist)/spkN;
    perc(4, ii) = length(northDist)/spkN;
end

[~, highWallIdx] = max(max(perc, [], 2));

end

%% 
% global nPosBins;
% distVec = linspace(0, 100, nPosBins);
% baseMat = repmat(distVec, nPosBins, 1);
% distMat(:, :, 1) = baseMat; % west
% distMat(:, :, 2) = rot90(baseMat, 1); % south
% distMat(:, :, 3) = rot90(baseMat, 2); % east
% distMat(:, :, 4) = rot90(baseMat, 3); % north
% bestScore = 0;
% 
% for ii = 1:4
%     
%     matrixSum = sum(nansum(ratemap.*distMat(:,:,ii)));
%     rateSum = sum(nansum(ratemap));
%     score = matrixSum/rateSum;
%     
%     if score >= bestScore
%         bestScore = score;
%         dir = ii;
%     end
% end