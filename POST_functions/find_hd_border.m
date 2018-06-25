%%
dbstop if error;
warning('off','all')
global m n count dt nPosBins nHdBins nSpeedBins lag boxSize;
m = 3; n = 4; % subplots
nPosBins = 50; nHdBins = 60; nSpeedBins = 10; lag = 500; dt = 0.04; boxSize = 100;

%% find_valid_trials
% files = subdir('C:\Users\elhan\Github\Cells\*.mat');
% saveIdx = 1;
% clearvars trialsToSave;
% 
% for ii = 1:length(files)
%     load(files(ii).name);
%     c = c;
%     vt = cVt;
%     
%     p = p;
%     s = p.S(c.session);
%     
%     % not NLG, not less than 10 minutes, and not less than 300 spikes
%     if strcmpi(p.nlgnlx, 'nlg') || ((s.end_time - s.start_time).*1e-6)/60 <= 10 || length(c.timestamps) <= 300
%         continue;
%     else
%         trialsToSave(saveIdx) = files(ii);
%         saveIdx = saveIdx + 1;
%     end
% end

%% find cells
% clearvars -except trialsToSave

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
    
    %% border
    [orgPosOccupancy, orgPosSpikes, ~, orgPosRatesSmooth, orgPosRatesSmoothUnclean] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx);
    
    [orgBorderScore, orgDir] = calculate_border_score(orgPosRatesSmooth);
    
    [orgBorderScoreUnclean, orgDirUnclean] = calculate_border_score(orgPosRatesSmoothUnclean);

    %% hd
    [orgHdOccupancy, ~, orgHdRates, ~] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);  
    orgHdScore = circ_r(hdBins', orgHdRates');
    shuffledHdScore = zeros(1, 1000);
    
    %% main functions
    POST_shuffle_hd_border;
    POST_plod_border;
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
