function [speedOccupancy, speedSpikes, speedRates] = calculate_speed_map(c, vt, cKeepIdx, vtKeepIdx)
global nSpeedBins dt;
speedBins = linspace(min(vt.speed(vtKeepIdx)), max(vt.speed(vtKeepIdx)), nSpeedBins);

speedOccupancy = histcounts(vt.speed(vtKeepIdx), speedBins);
speedSpikes = histcounts(c.speed(cKeepIdx), speedBins);

speedRates = speedSpikes./(speedOccupancy.*dt);
end