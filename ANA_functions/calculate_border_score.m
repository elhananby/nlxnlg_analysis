function [bestScore, dir] = calculate_border_score(ratemap)
global nPosBins;
distVec = linspace(0, 100, nPosBins);
baseMat = repmat(distVec, nPosBins, 1);
distMat(:, :, 1) = baseMat; % west
distMat(:, :, 2) = rot90(baseMat, 1); % south
distMat(:, :, 3) = rot90(baseMat, 2); % east
distMat(:, :, 4) = rot90(baseMat, 3); % north
bestScore = 0;

for ii = 1:4
    
    matrixSum = sum(nansum(ratemap.*distMat(:,:,ii)));
    rateSum = sum(nansum(ratemap));
    score = matrixSum/rateSum;
    
    if score >= bestScore
        bestScore = score;
        dir = ii;
    end
end

end