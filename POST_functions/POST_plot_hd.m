
hdThr = quantile(shuffledHdScore, 0.99);

if orgHdScore >= hdThr
    hdCells(hdIdx).cell = c.cell_number;
    hdCells(hdIdx).animal = c.animal;
    hdCells(hdIdx).score = orgHdScore;
    hdCells(hdIdx).rateMap = orgHdRates;
    hdIdx = hdIdx + 1;
    
    % save figure
    hdFig = figure('Units', 'normalized', 'Position', [0 0 0.3 0.75]);
    
    % smoothing
    win = hamming(10);
    win = win/sum(win);
    
    % hd over time
    subplot(3,2,[1 2])
    plotVtTimestamps = ((vt.timestamps-min(vt.timestamps)).*1e-6)/60;
    plotCTimestamps = ((c.timestamps - min(vt.timestamps)).*1e-6)/60;
    plot(plotVtTimestamps, vt.poshd, '-k', plotCTimestamps, c.poshd, '.r');
    title(sprintf('%i %i-%s %i\nTT%i C%i', c.cell_number, c.animal, c.animal_name, c.day, c.TT, c.cell_id));
    axis tight
    box off
    
    % occupancy
    occupancyHdScore = circ_r(hdBins', orgHdOccupancy');
    occupancyHdPval = circ_rtest(hdBins', orgHdOccupancy');
    occupancyHdMean = circ_mean(hdBins', orgHdOccupancy');
    
    subplot(3,2,3)
    polarplot(hdBins, cconv(orgHdOccupancy, win, nHdBins));
    title(sprintf('Occupancy\nr = %.2f p = %.2f mu = %0.f', occupancyHdScore, occupancyHdPval, rad2deg(wrapTo2Pi(occupancyHdMean))));
    hold on;
    polarplot([occupancyHdMean occupancyHdMean], [0 occupancyHdScore], 'r', 'LineWidth', 3);
    
    % rate
    hdScore = circ_r(hdBins', orgHdRates');
    hdPval = circ_rtest(hdBins', orgHdRates');
    hdMean = circ_mean(hdBins', orgHdRates');
    
    subplot(3,2,4)
    polarplot(hdBins, cconv(orgHdRates, win, nHdBins));
    title(sprintf('Rate\nr = %.2f p = %.2f mu = %.0f', hdScore, hdPval, rad2deg(wrapTo2Pi(hdMean))));
    
    % hd score
    hold on;
    polarplot([hdMean hdMean], [0 orgHdScore], 'r', 'LineWidth', 3);
    
    % shuffling
    subplot(3,2,[5 6])
    histogram(shuffledHdScore, 40);
    hold on;
    line([hdThr, hdThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    line([orgHdScore, orgHdScore], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    title(sprintf('Shuffled HD Scores\nthr = %.2f score = %.2f', hdThr, orgHdScore));
    xlabel('HD Score');
    ylabel('Count');
    
    % save figure
    hdFigFile = fullfile(figDir, sprintf('%i_%i_%i_hd', c.cell_number, c.animal, c.day));
    export_fig(hdFig, hdFigFile, '-jpg', '-native');
    close all;  
end