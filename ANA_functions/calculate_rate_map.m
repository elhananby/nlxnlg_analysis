function [posOccupancy, posSpikes, posRates, posRatesSmooth, varargout] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx)
global nPosBins boxSize dt;

binEdges = linspace(0, boxSize, nPosBins+1);

posOccupancy = histcounts2(vt.posx_c(vtKeepIdx), vt.posy_c(vtKeepIdx), binEdges, binEdges);

% find areas the animal spent less than % percent of the maximum time in
% and remove them
% timeOutliers = posOccupancy(:) < max(posOccupancy(:)) * 0.1;
% posOccupancy(timeOutliers) = NaN;

posSpikes = histcounts2(c.posx_c(cKeepIdx), c.posy_c(cKeepIdx), binEdges, binEdges);

posRates = posSpikes./(posOccupancy .* dt); % normalized to seconds to get frequency in Hz

rateOutliers = posRates >= nanmean(posRates(:))+6*nanstd(posRates(:));

posRates(isinf(posRates)) = NaN; % remove infs

varargout{1} = nanconv(posRates, fspecial('gaussian', 3*[3, 3], 3)); % original rate map

posRates(rateOutliers) = NaN; % clean outliers

posRatesSmooth = nanconv(posRates, fspecial('gaussian', 3*[3, 3], 3)); % smooth

end