function [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx)
global nPosBins boxSize dt;

binEdges = linspace(0, boxSize, nPosBins+1);

posOccupancy = histcounts2(vt.posx_c(vtKeepIdx), vt.posy_c(vtKeepIdx), binEdges, binEdges);

% find areas the animal spent less than % percent of the maximum time in
% and remove them
% timeOutliers = posOccupancy(:) < max(posOccupancy(:)) * 0.05;
% posOccupancy(timeOutliers) = NaN;

posSpikes = histcounts2(c.posx_c(cKeepIdx), c.posy_c(cKeepIdx), binEdges, binEdges);

posRates = posSpikes./(posOccupancy .* (dt*1e-6)); % normalized to seconds to get frequency in Hz

posRates(isinf(posRates)) = NaN;

posRatesSmooth = nanconv(posRates, fspecial('gaussian', 3*[1, 1], 1));

end