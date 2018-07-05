function [posOccupancy, posSpikes, posRates, posRatesSmooth, varargout] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx)
global nPosBins boxSize dt;

binEdges = linspace(0, boxSize, nPosBins+1);

posOccupancy = histcounts2(vt.posx_c(vtKeepIdx), vt.posy_c(vtKeepIdx), binEdges, binEdges);
posSpikes = histcounts2(c.posx_c(cKeepIdx), c.posy_c(cKeepIdx), binEdges, binEdges);

% outlier definition
posOccupancy(posOccupancy <= 1/dt) = 0; % remove all bins where the animal spent less than threshold

posRates = posSpikes ./ (posOccupancy.*dt); % calculate rate map normalized to the time in seconds

posRates(isinf(posRates)) = NaN; % remove infs

posRatesSmooth = nanconv(posRates, fspecial('gaussian', [3 3], 3), 'edge', 'nanout'); % smooth
end