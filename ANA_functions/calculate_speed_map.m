function [speedOccupancy, speedSpikes, speedRates, speedScore] = calculate_speed_map(c, vt)
global nSpeedBins dt;
speedBins = linspace(min(vt.speed), max(vt.speed), nSpeedBins);

speedOccupancy = histcounts(vt.speed, speedBins);
speedSpikes = histcounts(c.speed, speedBins);

speedRates = speedSpikes./(speedOccupancy.*dt);

speedScore = calculate_speed_score(c, vt);
end

function score = calculate_speed_score(c, vt)
global dt;
% correlation coefficient between instantaneous firing rate and speed
timebins = [vt.timestamps; (vt.timestamps(end) + dt)];
spikeTrain = histcounts(c.timestamps, timebins)';

% Smooths firing rate
filter = gaussmf(-4:4, [2 0]); filter = filter / sum(filter);
firingRate = spikeTrain / dt;
smoothFiringRate = conv(firingRate, filter, 'same');

score = corrcoef(smoothFiringRate, vt.speed); % calculate correlation
score = score(2); % get only the one we're interested in
end