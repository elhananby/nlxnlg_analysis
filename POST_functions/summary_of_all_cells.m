clearvars
close all
dbstop if error

files = subdir(fullfile('D:\', 'experiment_data', 'cells', '*.mat'));
files = cellstr(char(files.name));

cells = cell(length(files), 7);
% COLUMNS
%   1 - cell num, 2 - time, 3 - # spikes, 4 - lratio, 5 - isolationDist,
%   6 - behavior rayliehg, 7 - % ISI < 2ms

for ii_file = files'
    file = ii_file{1};
    
    disp(file);
    load(file);
    
    t = (C.S.pos.t(end) - C.S.pos.t(1))/60; % time in minutes
    ns = length(C.S.spk.t); % number of spikes
    lRatio = C.S.cluster_quality.Lratio; % L-Ratio
    isolationDist = C.S.cluster_quality.IsolationDist; % Isolation-distnace
    countRs = rayscore(C); % rayleigh score of BEHAVIOR (we want it to be as low as possible)
    ISI = diff(C.S.spk.Timestamps); % ISI
    perc_ISI_less_2ms = nnz(ISI .* 1e-3 <= 2)/length(ISI); % find ISI less than 2 ms
    
    [~, name, ~] = fileparts(file);
    number = regexp(name, '(\d+)_quail', 'tokens');
    number = str2double(number{1});
    
    idx = find(ismember(files, file));
    cells{idx, 1} = number;
    cells{idx, 2} = t;
    cells{idx, 3} = ns;
    cells{idx, 4} = lRatio;
    cells{idx, 5} = isolationDist;
    cells{idx, 6} = countRs;
    cells{idx, 7} = perc_ISI_less_2ms;
    
end

function rs = rayscore(C)
nHdBins = length(C.S.HD.time_phi);
hdBins = (0:2*pi/nHdBins:2*pi-2*pi/nHdBins)';

hdCounts = wrapTo2Pi(deg2rad(C.S.HD.time_phi))';
% convert coordinates from polar to rectangular form
x = hdCounts .* cos(hdBins);
y = hdCounts .* sin(hdBins);
sumLength = sum(sqrt(x.^2 + y.^2));

% computes sum vector
rx = sum(x);
ry = sum(y);
rLength = sqrt(rx^2 + ry^2);

rs = rLength / sumLength;
end