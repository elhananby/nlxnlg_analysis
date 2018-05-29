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