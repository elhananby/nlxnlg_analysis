function shiftedTimestamps = shift_timestamps(cTimestamps, vtTimestamps, varargin)


maxT = max(vtTimestamps);
minT = min(vtTimestamps);
maxJitter = maxT - minT;

switch nargin 
    case 2, jitter = maxJitter*rand;
    case 3, jitter = maxJitter*varargin{1};
end

shiftedTimestamps = mod(((cTimestamps-minT) + jitter), maxJitter) + minT;

end