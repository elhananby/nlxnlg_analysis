function ANA_basic_analysis(cell_list)
dbstop if error;
close all; clc;

global m n count dt nPosBins nHdBins nSpeedBins lag boxSize;
m = 3; n = 4; % subplots
nPosBins = 30; nHdBins = 60; nSpeedBins = 10; lag = 500;

%% enumerate over cells and load one cell at a time

for ii_cell = cell_list
    
    load(ii_cell{1});
    fprintf('Animal %i-%s\tDay %i\tExperiment %i\tSession %i\tTetrode %i\tCell %i\n',...
        p.animal, p.animal_name, p.day, p.experiment, c.session, c.TT, c.cell_id);

    % correct struct names for analysis
    c = c;
    vt = cVt;
    p = p;
    s = p.S(c.session);
    
    boxSize = p.arena_width_east_to_west;
    dt = mean(diff(vt.timestamps));

    fig = figure('Units', 'Normalized',...
        'Position', [0 0 1 1],...
        'Visible', 'off');
    
    %% cleaning thresholds
    vtKeepIdx = index_to_keep(vt, p, s);
    cKeepIdx = index_to_keep(c, p, s);
    
    %% calculate basic stuff
    % rate map
    [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx);
    
    % hd map
    [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);
    
    % speed map
    [speedOccupancy, speedSpikes, speedRates] = calculate_speed_score(c, vt, cKeepIdx, vtKeepIdx);
    
    % isi
    spikeISI = diff(c.timestamps.*1e-6); % in seconds
    
    % xcorr
    [rSpikeTrain, lagsSpikeTrain] = xcorr(c.spikeTrain, 500);
    
    % border score
    borderScore = calculate_border_score(posRates);
    
    %% plot all the basic stuff
    count = 1;
    
    %% position map
    subplot(m, n, count);
    plot(vt.posx_c, vt.posy_c, '-k', c.posx_c, c.posy_c, '.r');
    axis tight
    axis off
    
    count = count +1;
    
    %% rate map
    subplot(m, n, count);
    pcolor(posRatesSmooth');
    colormap(jet);
    shading flat;
    axis tight;
    axis off;
    colorbar;
    
    count = count + 1;
    
    %% spike shape
    subplot(m, n, count);
    
    % find channel with maximum value
    [~, maxChan] = max(max(mean(c.spikeShape, 3)));
    
    % plot all spike lines
    plot(squeeze(c.spikeShape(:, maxChan, :)), 'Color', [220 220 220 50]./255, 'LineWidth', .5);
    hold on;
    
    % plot only mean of spike
    plot(mean(squeeze(c.spikeShape(:, maxChan, :)), 2), 'k', 'LineWidth', 1);
    
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    box off;
    
    count = count + 1;
    
    %% ISI
    subplot(m, n, count);
    x = logspace(-4, 1, 100);
    histogram(spikeISI, x);
    set(gca, 'xscale', 'log');
    hold on
    line([0.002 0.002], [0 max(histcounts(log(spikeISI), 100))], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
    xlabel('Time (ms)');
    ylabel('ISI Count');
    title(sprintf('%i# Spikes\n%.2fHz',...
        length(c.timestamps),...
        length(c.timestamps)/((vt.timestamps(end)-vt.timestamps(1))*1e-6)));
    box off
    count = count + 1;
    
    %% head direction time plot
    subplot(m, n, count);
    angleDiff = diff(vt.poshd);
    correctedPosHd = vt.poshd;
    correctedPosHd(abs(angleDiff) >= 2*pi/1.8) = NaN;
    plot(vt.timestamps, correctedPosHd, '-k', c.timestamps, c.poshd, '.r');
    count = count + 1;
    
    %% head direction polar plot
    % smoothing
    win = hamming(10);
    win = win/sum(win);
    
    % polar plot
    subplot(m, n, count);
    hdBins = linspace(-pi, pi, nHdBins);
    
    % occupancy
    polarplot(hdBins, cconv(hdOccupancy/max(hdOccupancy), win, nHdBins));
    
    hold on;
    % rate
    polarplot(hdBins, cconv(hdRates/max(hdRates), win, nHdBins));
    
    % basic calculations
    hdMean = circ_mean(hdBins', hdRates');
    hdScoreP = circ_rtest(hdBins', hdRates');
    hdScore = circ_r(hdBins', hdRates');
    
    % hd score
    polarplot([hdMean hdMean], [0 hdScore], 'r', 'LineWidth', 3);
    
    title(sprintf('Rayleigh = %.2f \t p = %d \n Mean Direction = %1.2f',...
        hdScore, hdScoreP,rad2deg(hdMean)), 'FontSize', 8);
    count = count + 1;
    
    % histogram
    subplot(m, n, count);
    hdBins = -pi : 2*pi/nHdBins : pi - pi/nHdBins;
    yyaxis left;
    area(hdBins, cconv(hdOccupancy, win, nHdBins), 'FaceAlpha', .5, 'LineStyle', 'none');
    ylabel('Time (sec)');
    
    yyaxis right;
    area(hdBins, cconv(hdRates, win, nHdBins), 'FaceAlpha', .5, 'LineStyle', 'none');
    xlabel('Head Direction (rad)');
    ylabel('Firing Rate (Hz)');
    
    axis tight
    box off
    
    hdOccupancyP = circ_rtest(hdBins', hdOccupancy');
    hdOccupancyScore = circ_r(hdBins', hdOccupancy');
    
    title(sprintf('Normalized Rayleigh Score: \n %.2f Rate / %.2f Behavior = %.2f',...
        hdScore, hdOccupancyScore, hdScore/hdOccupancyScore));
    count = count + 1;
    
    % head direction shuffling
    subplot(m, n, count);
    shuffledHdScore = zeros(1, 1000);
    for ii_shuffle = 1:1000
        
        % shuffling and calculations
        shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
        iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
        iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
        
        % variable of interest
        [~, ~, shuffledHdRates, ~] = calculate_hd_map(iC, vt, iCKeepIdx, vtKeepIdx);
        hdBins = linspace(-pi, pi, nHdBins);
        shuffledHdScore(ii_shuffle) = circ_r(hdBins', shuffledHdRates');
    end
    
    hdThr = quantile(shuffledHdScore, 0.99);
    histogram(shuffledHdScore, 40);
    hold on;
    line([hdThr, hdThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    hdScore = circ_r(hdBins', hdRates');
    line([hdScore hdScore], ylim, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
    
    title('Shuffled HD Scores');
    xlabel('HD Score');
    ylabel('Count');
    
    count = count + 1;
    
    %% shuffle border score
    subplot(m, n, count);
    shuffledBorderScore = zeros(1, 1000);
    
    for ii_shuffle = 1:1000
        
        % shuffling and calculations
        shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
        iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
        iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
        
        % calculate variable
        [~, ~, ~, posRatesSmooth] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx); % calculate rate map
        shuffledBorderScore(ii_shuffle) = calculate_border_score(posRatesSmooth); % calculate border score
    end
    
    borderThr = quantile(shuffledBorderScore, 0.99);
    histogram(shuffledBorderScore, 40);
    hold on;
    line([borderThr, borderThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([borderScore, borderScore], ylim, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
    
    title('Shuffled Border Scores');
    xlabel('Border Score');
    ylabel('Count');
    
    count = count + 1;
    
    %% speed
    subplot(m ,n, count);
    speedBins = linspace(min(vt.speed(vtKeepIdx)), max(vt.speed(vtKeepIdx)), nSpeedBins-1);
    yyaxis left;
    area(speedBins, speedOccupancy, 'FaceAlpha', .5, 'LineStyle', 'none');
    ylabel('Time (sec)');
    
    yyaxis right;
    area(speedBins, speedRates, 'FaceAlpha', .5, 'LineStyle', 'none');
    xlabel('Speed (cm/s)');
    ylabel('Firing Rate (Hz)');
    
    axis tight
    box off
    count = count + 1;
    
    %% title
    suptitle(sprintf('%s \t Animal %i \t Day %i \n Experiment %i \t Session %i - %s \n TT %i \t Cell %i',...
        p.nlgnlx,...
        p.animal,...
        p.day,...
        p.experiment,...
        c.session,...
        p.S(c.session).session,...
        c.TT,...
        c.cell_id));
    
    % save
    [filepath_fig, ~, ~] = fileparts(ii_cell{1});
    filename_fig = sprintf('%s_Animal%i_Day%i_Exp%i_Session%i_TT%i_C%i',...
        p.nlgnlx,...
        p.animal,...
        p.day,...
        p.experiment,...
        c.session,...
        c.TT,...
        c.cell_id);
    
    savefig(fig, fullfile(filepath_fig, filename_fig), 'compact');
    saveas(fig, fullfile(filepath_fig, [filename_fig '.png']), 'png');
    close all;
    fprintf('\t\tSaved %s\n', fullfile(filepath_fig, filename_fig));
end % cells

end % functions


%% calculation functions
function [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx)
global nPosBins boxSize dt;

binEdges = linspace(0, boxSize, nPosBins+1);

posOccupancy = histcounts2(vt.posx_c(vtKeepIdx), vt.posy_c(vtKeepIdx), binEdges, binEdges);

% find areas the animal spent less than % percent of the maximum time in
% and remove them
% timeOutliers = posOccupancy(:) < max(posOccupancy(:)) * 0.05;
% posOccupancy(timeOutliers) = NaN;

posSpikes = histcounts2(c.posx_c(cKeepIdx), c.posy_c(cKeepIdx), binEdges, binEdges);

posRates = posSpikes./(posOccupancy .* (dt*1e-6)); % normalized to seconds to get frequency in Hz

posRates(isinf(posRates)) = NaN;

posRatesSmooth = nanconv(posRates, fspecial('gaussian', 3*[1, 1], 1));

end

function [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx)
global nHdBins dt;
hdBins = -pi : 2*pi/nHdBins : pi;

hdOccupancy = histcounts(vt.poshd(vtKeepIdx), hdBins);
hdSpikes = histcounts(c.poshd(cKeepIdx), hdBins);
hdRates = hdSpikes ./ (hdOccupancy .* (dt * 1e-6));

% convert coordinates from polar to rectangular form
hdBins = -pi : 2*pi/nHdBins : pi - pi/nHdBins;
x = hdRates .* cos(hdBins);
y = hdRates .* sin(hdBins);
sumLength = sum(sqrt(x.^2 + y.^2));

% computes sum vector
rx = sum(x);
ry = sum(y);
rLength = sqrt(rx^2 + ry^2);

hdScore = rLength / sumLength;

end

function [speedOccupancy, speedSpikes, speedRates] = calculate_speed_score(c, vt, cKeepIdx, vtKeepIdx)
global nSpeedBins dt;
speedBins = linspace(min(vt.speed(vtKeepIdx)), max(vt.speed(vtKeepIdx)), nSpeedBins);

speedOccupancy = histcounts(vt.speed(vtKeepIdx), speedBins);
speedSpikes = histcounts(c.speed(cKeepIdx), speedBins);

speedRates = speedSpikes./(speedOccupancy.*(dt*1e-6));
end

function borderScore = calculate_border_score(rateMap)
maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);

% Calculates CM
nFields = max(modifiedRateMap(:));
maxDistance = 0;
longestFieldId = 0;
maxIdx = 0;
for i = 1:nFields
    % Determines number of adjacent pixels to a wall in a firing field
    [fieldRow, fieldCol] = find(modifiedRateMap == i);
    
    nFieldPixels = length(fieldRow);
    % Fields defined as needing to have more than 6 pixels
    if (nFieldPixels <= 6)
        %         modifiedRateMap(modifiedRateMap == i) = 0;
        continue;
    end
    
    % Calculates length of pixels along each wall
    rowLeftPixels = fieldRow(fieldCol == 1); % pixels along vertical left wall
    leftDistance = abs(max(rowLeftPixels) - min(rowLeftPixels));
    if isempty(leftDistance)
        leftDistance = 0;
    end
    
    rowRightPixels = fieldRow(fieldCol == length(rateMap)); % pixels along vertical right wall
    rightDistance = abs(max(rowRightPixels) - min(rowRightPixels));
    if isempty(rightDistance)
        rightDistance = 0;
    end
    
    colTopPixels = fieldCol(fieldRow == 1); % pixels along horizontal top wall
    topDistance = abs(max(colTopPixels) - min(colTopPixels));
    if isempty(topDistance)
        topDistance = 0;
    end
    
    colBottomPixels = fieldCol(fieldRow == length(rateMap)); % pixels along horizontal bottom wall
    bottomDistance = abs(max(colBottomPixels) - min(colBottomPixels));
    if isempty(bottomDistance)
        bottomDistance = 0;
    end
    
    % Determines maximal length of wall touching the firing field
    [fieldMaxDistance, idx] = max([leftDistance rightDistance topDistance bottomDistance]);
    
    if fieldMaxDistance > maxDistance
        longestFieldId = i;
        maxIdx = idx;
        maxDistance = fieldMaxDistance;
    end
end

CM = maxDistance;

% Calculates DM
% Commented out lines used for Method 2
% [row, col] = find(modifiedRateMap == longestField);
[row, col] = find(modifiedRateMap ~= 0);
nPixels = length(row);
totalDistance = 0;
for i = 1:nPixels
    %     if maxIdx == 1 || maxIdx == 3
    %         wall = 1;
    %     else
    %         wall = length(rateMap);
    %     end
    %
    %     if maxIdx == 1 || maxIdx == 2
    %         distance = abs(col(i) - wall);
    %     else
    %         distance = abs(row(i) - wall);
    %     end
    
    distance = min([abs(col(i) - 1) abs(col(i) - length(rateMap)) ...
        abs(row(i) - 1) abs(row(i) - length(rateMap))]);
    
    totalDistance = totalDistance + distance;
end

DM = totalDistance / nPixels;
if (isnan(DM))
    DM = 0;
end

borderScore = (CM - DM)/(CM + DM);

end

%% support functions
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

function iC = interpolate_shifted_values(ts, vt)

iC.timestamps = ts;
iC.posx       = interp1(vt.timestamps, vt.posx, ts);
iC.posx2      = interp1(vt.timestamps, vt.posx2, ts);
iC.posy       = interp1(vt.timestamps, vt.posy, ts);
iC.posy2      = interp1(vt.timestamps, vt.posy2, ts);
iC.posx_c     = interp1(vt.timestamps, vt.posx_c, ts);
iC.posy_c     = interp1(vt.timestamps, vt.posy_c, ts);
iC.poshd      = interp1(vt.timestamps, vt.poshd, ts);
iC.vx         = interp1(vt.timestamps, vt.vx, ts);
iC.vy         = interp1(vt.timestamps, vt.vy, ts);
iC.speed      = interp1(vt.timestamps, vt.speed, ts);

end

function keepIdx = index_to_keep(arr, p, s)

session_time = arr.timestamps <= s.start_time & arr.timestamps >= s.end_time; % get only session times
speed_threshold = arr.speed >= 50 | arr.speed <= 1; % get speed threshold
no_behavior = PRE_throw_away_times(arr.timestamps, p.throw_away_times); % remove throw away times
keepIdx = ~(session_time | speed_threshold | no_behavior); % NOT throw away to get indices to keep

end