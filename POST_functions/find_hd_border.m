%%
clearvars
dbstop if error;
warning('off','all')
global m n count dt nPosBins nHdBins nSpeedBins lag boxSize;
m = 3; n = 4; count = 1; % subplots
nPosBins = 50; nHdBins = 60; nSpeedBins = 10; lag = 500; dt = 0.04; boxSize = 100;

%% find_valid_trials
files = subdir('C:\Users\elhan\Github\Cells\*.mat');
saveIdx = 1;
clearvars trialsToSave;
lineCount = 0;
if ~exist('trialsToSave.mat', 'file')
    for ii = 1:length(files)
        load(files(ii).name);
        c = c;
        vt = cVt;
        
        p = p;
        s = p.S(c.session);
        
        % not NLG, not less than 10 minutes, and not less than 300 spikes
        if strcmpi(p.nlgnlx, 'nlg') || ((s.end_time - s.start_time).*1e-6)/60 <= 10 || length(c.timestamps) <= 300
            continue;
        else
            trialsToSave(saveIdx) = files(ii);
            saveIdx = saveIdx + 1;
        end
        fprintf(repmat('\b', 1, lineCount));
        lineCount = fprintf('%i/%i', ii, length(files));
    end
    save('trialsToSave.mat', 'trialsToSave');
else, load('trialsToSave.mat');
end
%% find cells
figDir = fullfile('C:\Users\elhan\Github\Cells');

lineCount = 0;
totTime = tic;

borderCells = struct;
hdCells = struct;

borderIdx = 1;
hdIdx = 1;

hdBins = linspace(-pi, pi, nHdBins);
    
for ii = 1:length(trialsToSave)

    innTime = tic;
    
    %% Load cell data
    load(trialsToSave(ii).name);
    c = c;
    vt = cVt;
    p = p;
    s = p.S(c.session);
    
    vtKeepIdx = index_to_keep(vt, p, s);
    cKeepIdx = index_to_keep(c, p, s);
    
%     %% border
%     [orgPosOccupancy, orgPosSpikes, ~, orgPosRatesSmooth, orgPosRatesSmoothUnclean] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx);
%     
%     [orgBorderScore, orgDir] = calculate_border_score(orgPosRatesSmooth);
%     
%     [orgBorderScoreUnclean, orgDirUnclean] = calculate_border_score(orgPosRatesSmoothUnclean);

    %% hd
    % total
    [orgHdOccupancy, ~, orgHdRates, ~] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);  
    orgHdScore = circ_r(hdBins', orgHdRates');
    
    % first half
    halfToKeep{1} = find(vt.timestamps <= vt.timestamps(ceil(length(vt.timestamps)/2)))';
    halfToKeep{2} = find(c.timestamps <= vt.timestamps(ceil(length(vt.timestamps)/2)))';
    
    [hdOccupancy1, ~, hdRates1, ~] = calculate_hd_map(c, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    orgHdScore1 = circ_r(hdBins', hdRates1');

    % second half
    halfToKeep{1} = find(vt.timestamps > vt.timestamps(ceil(length(vt.timestamps)/2)))';
    halfToKeep{2} = find(c.timestamps > vt.timestamps(ceil(length(vt.timestamps)/2)))';
    
    [hdOccupancy2, ~, hdRates2, ~] = calculate_hd_map(c, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    orgHdScore2 = circ_r(hdBins', hdRates2');
    
    % slow
    halfToKeep{1} = find(vt.speed <= 1)';
    halfToKeep{2} = find(c.speed <= 1)';
    
    [hdOccupancySlow, ~, hdRatesSlow, ~] = calculate_hd_map(c, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    orgHdScoreSlow = circ_r(hdBins', hdRatesSlow');
    
    % fast
    halfToKeep{1} = find(vt.speed > 1)';
    halfToKeep{2} = find(c.speed > 1)';
    
    [hdOccupancyFast, ~, hdRatesFast, ~] = calculate_hd_map(c, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    orgHdScoreFast = circ_r(hdBins', hdRatesFast');
    
    shuffledHdScore = zeros(1, 1000);
    shuffledHdScore1 = zeros(1, 1000);
    shuffledHdScore2 = zeros(1, 1000);
    shuffledHdScoreSlow = zeros(1, 1000);
    shuffledHdScoreFast = zeros(1, 1000);
    
    %% main functions
    POST_shuffle_hd_border;
%     POST_plod_border;
    POST_plot_hd;
    
    %% progress bar
    fprintf(repmat('\b', 1, lineCount));
    
    lineCount = fprintf('%i/%i %.2f %.2f%% %.2f Seconds',...
        ii,...
        length(trialsToSave),...
        toc(innTime),...
        ii*100/length(trialsToSave),...
        toc(totTime));

end

save('hdcells.mat', 'hdCells');
save('bordercells.mat', 'borderCells');

fprintf('\nDone');

%% plot heatmap of rates
heatRateMap = figure();
hdBins = linspace(-pi, pi, 60);
heatRateMap = [];
idx = 1;

for ii = 1:length(hdCells)
    rateMap = hdCells(ii).rateMap;
    if hdCells(ii).score >= 0.2
        mu(idx) = circ_mean(hdBins', rateMap');
        heatRateMap(idx, :) = rateMap./max(rateMap);
        idx = idx + 1;
    end
end
[B, I] = sort(mu);
heatRateMapSorted = heatRateMap(I, :);
imagesc(heatRateMapSorted);
xlim([0 60]);
xticklabels(0:60:360);
hdRateFigFile = fullfile(figDir, 'hdHeatMap');
export_fig(heatRateMap, hdRateFigFile, '-jpg');