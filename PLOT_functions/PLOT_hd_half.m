subplot(m, n, count);

hdBinsCalculate = -pi : 2*pi/nHdBins : pi;
hdBinsPlot = linspace(-pi, pi, 60);

midTime = vt.timestamps(ceil(length(vt.timestamps)/2));

hdRates1 = histcounts(c.poshd(c.timestamps <= midTime), hdBinsCalculate) ./...
    (histcounts(vt.poshd(vt.timestamps <= midTime), hdBinsCalculate) .* dt );
hdMean1 = circ_mean(hdBinsPlot', hdRates1');
hdScore1 = circ_r(hdBinsPlot', hdRates1');

hdRates1(isnan(hdRates1)) = 0;

hdRates2 = histcounts(c.poshd(c.timestamps > midTime), hdBinsCalculate) ./...
    (histcounts(vt.poshd(vt.timestamps > midTime), hdBinsCalculate) .* dt);
hdMean2 = circ_mean(hdBinsPlot', hdRates2');
hdScore2 = circ_r(hdBinsPlot', hdRates2');

hdRates2(isnan(hdRates2)) = 0;

p1 = polarplot(hdBinsPlot, cconv(hdRates1./max(hdRates1), win, nHdBins), 'Color', [0 0.4470 0.7410]);
hold on
polarplot([hdMean1 hdMean1], [0 hdScore1], 'Color', [0 0.4470 0.7410], 'LineWidth', 3); 
p2 = polarplot(hdBinsPlot, cconv(hdRates2./max(hdRates2), win, nHdBins), 'Color', [0.8500 0.3250 0.0980]);
polarplot([hdMean2 hdMean2], [0 hdScore2], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 3);

legend([p1 p2], {'First', 'Second'}, 'location', 'best');
title('Halves Analysis');

count = count + 1;
