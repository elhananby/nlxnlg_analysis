function [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(varargin)
%% CALCULATE_RATE_MAP used to claculate the behavior rate map.
%  INPUT:   c - cell structure; vt = video structre OR
%           vt.posx, vt.posy, c.posx, c.posy - (IN THAT ORDER) location
%           data
global nPosBins boxSize dt;

%% input processing
if nargin == 2
    c = varargin{1};
    vt = varargin{2};
elseif nargin == 4
    vt.posx_c = varargin{1};
    vt.posy_c = varargin{2};
    c.posx_c = varargin{3};
    c.posy_c = varargin{4};
else
    error('Error. Wrong number of inputs; either 2 or 4');
end

%% main
binEdges = 0:boxSize/nPosBins:boxSize;

posOccupancy = histcounts2(vt.posx_c, vt.posy_c, binEdges, binEdges);
posSpikes = histcounts2(c.posx_c, c.posy_c, binEdges, binEdges);


posRates = posSpikes ./ (posOccupancy.*dt); % calculate rate map normalized to the time in seconds

posRates(isinf(posRates)) = NaN; % remove infs

% outlier definition (mean + 5 std)
posRates(posRates >= nanmean(posRates(:)) + 5*nanstd(posRates(:))) = NaN;

posRatesSmooth = nanconv(posRates, fspecial('gaussian', 3*[3 3], 3)); % smooth
% end