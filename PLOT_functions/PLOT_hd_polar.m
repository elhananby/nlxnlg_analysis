%% head direction polar plot

% polar plot
subplot(m, n, count);

% smoothing
win = hamming(10);
win = win/sum(win);

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
hdHaScore = circ_otest(hdBins', [], hdRates');

% hd score
polarplot([hdMean hdMean], [0 hdScore], 'r', 'LineWidth', 3);

title(sprintf('Rayleigh = %.2f \t p = %d \n Mean Direction = %1.2f',...
    hdScore, hdScoreP, rad2deg(hdMean)), 'FontSize', 8);

count = count + 1;