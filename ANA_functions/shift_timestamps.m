function shiftedTimestamps = shift_timestamps(cTimestamps, vtTimestamps)

dt = mean(diff(vtTimestamps));
timebins = [vtTimestamps; (vtTimestamps(end) + dt)];
spikeTrain = histcounts(cTimestamps, timebins)';

minShift = ceil(20/dt); % min shift is 20 s
maxShift = length(spikeTrain)-(20/dt); % max shift is length of trial minus 20 s
randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by

shiftedSpiketrain = circshift(spikeTrain, randShift); % shifted spiketrain

shiftedTimestamps = RunLength(vtTimestamps, shiftedSpiketrain);

end