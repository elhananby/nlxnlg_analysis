dbstop if error;
m = 3;
n = 5;
hdThr = quantile(shuffledHdScore, 0.99);
hdThr1 = quantile(shuffledHdScore1, .99);
hdThr2 = quantile(shuffledHdScore2, .99);
hdThrSlow = quantile(shuffledHdScoreSlow, .99);
hdThrFast = quantile(shuffledHdScoreFast, .99);

if orgHdScore >= hdThr && orgHdScore >= .2
    hdCells(hdIdx).cell = c.cell_number;
    hdCells(hdIdx).animal = c.animal;
    hdCells(hdIdx).score = orgHdScore;
    hdCells(hdIdx).rateMap = orgHdRates;
    hdIdx = hdIdx + 1;
    
    % save figure
    hdFig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
    % smoothing
    win = hamming(10);
    win = win/sum(win);
    
    %% total
    subplot(m,n,1)
    plot(vt.posx_c, vt.posy_c, '-k', c.posx_c, c.posy_c, '.r');
    axis tight
    axis off
    
    subplot(m, n, 2)
    plot(mean(c.spikeShape(:, c.TT, :), 3))
    axis tight
    
    % hd over time
    subplot(m,n,[3 5])
    plotVtTimestamps = ((vt.timestamps-min(vt.timestamps)).*1e-6)/60;
    plotCTimestamps = ((c.timestamps - min(vt.timestamps)).*1e-6)/60;
    plot(plotVtTimestamps, rad2deg(wrapTo2Pi(vt.poshd)), '-k', plotCTimestamps, rad2deg(wrapTo2Pi(c.poshd)), '.r');
    title(sprintf('%i %i-%s %i\nTT%i C%i', c.cell_number, c.animal, c.animal_name, c.day, c.TT, c.cell_id));
    axis tight
    box off
    
    % occupancy
    occupancyHdScore = circ_r(hdBins', orgHdOccupancy');
    occupancyHdPval = circ_rtest(hdBins', orgHdOccupancy');
    occupancyHdMean = circ_mean(hdBins', orgHdOccupancy');
    
    subplot(m,n,6)
    polarplot(hdBins, cconv(orgHdOccupancy./max(orgHdOccupancy), win, nHdBins));
%     title(sprintf('Occupancy\nr = %.2f p = %.3f mu = %0.f', occupancyHdScore, occupancyHdPval, rad2deg(wrapTo2Pi(occupancyHdMean))));
    hold on;
    % rate
    hdScore = circ_r(hdBins', orgHdRates');
    hdPval = circ_rtest(hdBins', orgHdRates');
    hdMean = circ_mean(hdBins', orgHdRates');
    
    polarplot(hdBins, cconv(orgHdRates./max(orgHdRates), win, nHdBins));
    title(sprintf('Total n = %i\nr = %.2f p = %.3f mu = %.0f',...
        length(c.timestamps), hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    % hd score
    polarplot([hdMean hdMean], [0 orgHdScore], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(m, n, 11)
    histogram(shuffledHdScore, 40);
    hold on;
    line([hdThr, hdThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScore, orgHdScore], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThr, orgHdScore));
    xlabel('HD Score');
    ylabel('Count');
    
    orgvt = vt;
    orgc = c;
     
    %% first half
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    halfToKeep{1} = find(orgvt.timestamps <= orgvt.timestamps(ceil(length(orgvt.timestamps)/2)))';
    halfToKeep{2} = find(orgc.timestamps <= orgvt.timestamps(ceil(length(orgvt.timestamps)/2)))';
    vt = structfun(@(x) x(halfToKeep{1}), orgvt, 'UniformOutput', false);
    c = rmfield(orgc, {'cell_number'; 'animal'; 'animal_name'; 'day'; 'experiment'; 'session'; 'TT'; 'cell_id';...
        'spikeTrain'; 'firingRate'; 'smoothFiringRate'; 'spikeShape'});
    c = structfun(@(x) x(halfToKeep{2}), c, 'UniformOutput', false); 
        
    % occupancy
    subplot(m,n,7)
    occupancyHdScore = circ_r(hdBins', hdOccupancy1');
    occupancyHdPval = circ_rtest(hdBins', hdOccupancy1');
    occupancyHdMean = circ_mean(hdBins', hdOccupancy1');
    
    polarplot(hdBins, cconv(hdOccupancy1./max(hdOccupancy1), win, nHdBins));
    hold on;

    % rate
    hdScore = circ_r(hdBins', hdRates1');
    hdPval = circ_rtest(hdBins', hdRates1');
    hdMean = circ_mean(hdBins', hdRates1');
    
    polarplot(hdBins, cconv(hdRates1./max(hdRates1), win, nHdBins));
    title(sprintf('First Half n = %i\nr = %.2f p = %.3f mu = %.0f', length(halfToKeep{2}), hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    polarplot([hdMean hdMean], [0 orgHdScore1], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(m, n, 12)
    histogram(shuffledHdScore1, 40);
    hold on;
    line([hdThr1, hdThr1], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScore1, orgHdScore1], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThr1, orgHdScore1));
    xlabel('HD Score');
    ylabel('Count');
    
    
    %% second half
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    halfToKeep{1} = find(orgvt.timestamps > orgvt.timestamps(ceil(length(orgvt.timestamps)/2)))';
    halfToKeep{2} = find(orgc.timestamps > orgvt.timestamps(ceil(length(orgvt.timestamps)/2)))';
    
    vt = structfun(@(x) x(halfToKeep{1}), orgvt, 'UniformOutput', false);
    c = rmfield(orgc, {'cell_number'; 'animal'; 'animal_name'; 'day'; 'experiment'; 'session'; 'TT'; 'cell_id';...
        'spikeTrain'; 'firingRate'; 'smoothFiringRate'; 'spikeShape'});
    c = structfun(@(x) x(halfToKeep{2}), c, 'UniformOutput', false); 

    subplot(m,n,8)
    % occupancy
    occupancyHdScore = circ_r(hdBins', hdOccupancy2');
    occupancyHdPval = circ_rtest(hdBins', hdOccupancy2');
    occupancyHdMean = circ_mean(hdBins', hdOccupancy2');

    polarplot(hdBins, cconv(hdOccupancy2./max(hdOccupancy2), win, nHdBins));
    hold on;

    % rate
    hdScore = circ_r(hdBins', hdRates2');
    hdPval = circ_rtest(hdBins', hdRates2');
    hdMean = circ_mean(hdBins', hdRates2');
    
    polarplot(hdBins, cconv(hdRates2./max(hdRates2), win, nHdBins));
    title(sprintf('Second Half n = %i\nr = %.2f p = %.3f mu = %.0f', length(halfToKeep{2}), hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    polarplot([hdMean hdMean], [0 orgHdScore2], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(m, n, 13)
    histogram(shuffledHdScore2, 40);
    hold on;
    line([hdThr2, hdThr2], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScore2, orgHdScore2], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThr2, orgHdScore2));
    xlabel('HD Score');
    ylabel('Count');
    
    %% slow (<=1cm/s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    halfToKeep{1} = find(orgvt.speed <= 1)';
    halfToKeep{2} = find(orgc.speed <= 1)';
    
    vt = structfun(@(x) x(halfToKeep{1}), orgvt, 'UniformOutput', false);
    c = rmfield(orgc, {'cell_number'; 'animal'; 'animal_name'; 'day'; 'experiment'; 'session'; 'TT'; 'cell_id';...
        'spikeTrain'; 'firingRate'; 'smoothFiringRate'; 'spikeShape'});
    c = structfun(@(x) x(halfToKeep{2}), c, 'UniformOutput', false); 
    
    % occupancy
    subplot(m,n,9)
    occupancyHdScore = circ_r(hdBins', hdOccupancySlow');
    occupancyHdPval = circ_rtest(hdBins', hdOccupancySlow');
    occupancyHdMean = circ_mean(hdBins', hdOccupancySlow');
   
    polarplot(hdBins, cconv(hdOccupancySlow./max(hdOccupancySlow), win, nHdBins));
    hold on;

    % rate
    hdScore = circ_r(hdBins', hdRatesSlow');
    hdPval = circ_rtest(hdBins', hdRatesSlow');
    hdMean = circ_mean(hdBins', hdRatesSlow');
    
    polarplot(hdBins, cconv(hdRatesSlow./max(hdRatesSlow), win, nHdBins));
    title(sprintf('Slow (<= 1cm/s) n = %i\nr = %.2f p = %.3f mu = %.0f', length(halfToKeep{2}), hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    polarplot([hdMean hdMean], [0 orgHdScoreSlow], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(m, n, 14)
    histogram(shuffledHdScoreSlow, 40);
    hold on;
    line([hdThrSlow, hdThrSlow], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScoreSlow, orgHdScoreSlow], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThrSlow, orgHdScoreSlow));
    xlabel('HD Score');
    ylabel('Count');
    
    %% fast (>1cm/s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    halfToKeep{1} = find(orgvt.speed > 1)';
    halfToKeep{2} = find(orgc.speed > 1)';
    
    vt = structfun(@(x) x(halfToKeep{1}), orgvt, 'UniformOutput', false);
    c = rmfield(orgc, {'cell_number'; 'animal'; 'animal_name'; 'day'; 'experiment'; 'session'; 'TT'; 'cell_id';...
        'spikeTrain'; 'firingRate'; 'smoothFiringRate'; 'spikeShape'});
    c = structfun(@(x) x(halfToKeep{2}), c, 'UniformOutput', false); 
    
        subplot(m,n,10)
    % occupancy
    occupancyHdScore = circ_r(hdBins', hdOccupancyFast');
    occupancyHdPval = circ_rtest(hdBins', hdOccupancyFast');
    occupancyHdMean = circ_mean(hdBins', hdOccupancyFast');
    
    polarplot(hdBins, cconv(hdOccupancyFast./max(hdOccupancyFast), win, nHdBins));
    hold on;

    % rate
    hdScore = circ_r(hdBins', hdRatesFast');
    hdPval = circ_rtest(hdBins', hdRatesFast');
    hdMean = circ_mean(hdBins', hdRatesFast');
    
    polarplot(hdBins, cconv(hdRatesFast./max(hdRatesFast), win, nHdBins));
    title(sprintf('Fast (>1cm/s) n = %i\nr = %.2f p = %.3f mu = %.0f', length(halfToKeep{2}), hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    polarplot([hdMean hdMean], [0 orgHdScoreFast], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(m, n, 15)
    histogram(shuffledHdScore1, 40);
    hold on;
    line([hdThr1, hdThr1], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScore1, orgHdScore1], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThr1, orgHdScore1));
    xlabel('HD Score');
    ylabel('Count');
    
    %% save figure
    hdFigFile = fullfile(figDir, sprintf('%i_%i_%i_hd', orgc.cell_number, orgc.animal, orgc.day));
    export_fig(hdFig, hdFigFile, '-jpg', '-native');
    savefig(hdFig, hdFigFile);
    close all;  
end