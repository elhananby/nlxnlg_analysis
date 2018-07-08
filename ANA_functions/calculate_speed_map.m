function [speedOccupancy, speedSpikes, speedRates] = calculate_speed_map(c, vt)
global nSpeedBins dt;
speedBins = linspace(min(vt.speed), max(vt.speed), nSpeedBins);

speedOccupancy = histcounts(vt.speed, speedBins);
speedSpikes = histcounts(c.speed, speedBins);

speedRates = speedSpikes./(speedOccupancy.*dt);
end