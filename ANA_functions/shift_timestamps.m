function shiftedTimestamps = shift_timestamps(cTimestamps, vtTimestamps)

maxT = max(vtTimestamps);
minT = min(vtTimestamps);
maxJitter = maxT - minT;

jitter = maxJitter*rand;
shiftedTimestamps = mod(((cTimestamps-minT) + jitter), maxJitter) + minT;

end