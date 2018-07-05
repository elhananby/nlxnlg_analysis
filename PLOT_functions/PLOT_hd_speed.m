subplot(m, n, count);

hdBinsCalculate = -pi : 2*pi/nHdBins : pi;
hdBinsPlot = linspace(-pi, pi, 60);

vt.speed(vt.speed>50) = 50;

medSpeed = median(vt.speed);

hdRatesSlow = histcounts(c.poshd(c.speed <= medSpeed), hdBinsCalculate) ./...
    (histcounts(vt.poshd(vt.speed <= medSpeed), hdBinsCalculate) .* dt );
hdMeanSlow = circ_mean(hdBinsPlot', hdRatesSlow');
hdScoreSlow = circ_r(hdBinsPlot', hdRatesSlow');

hdRatesSlow(isnan(hdRatesSlow)) = 0;

hdRatesFast = histcounts(c.poshd(c.speed > medSpeed), hdBinsCalculate) ./...
    (histcounts(vt.poshd(vt.speed > medSpeed), hdBinsCalculate) .* dt);
hdMeanFast = circ_mean(hdBinsPlot', hdRatesFast');
hdScoreFast = circ_r(hdBinsPlot', hdRatesFast');

hdRatesFast(isnan(hdRatesFast)) = 0;

p1 = polarplot(hdBinsPlot, cconv(hdRatesSlow./max(hdRatesSlow), win, nHdBins), 'Color', [0 0.4470 0.7410]);
hold on
polarplot([hdMeanSlow hdMeanSlow], [0 hdScoreSlow], 'Color', [0 0.4470 0.7410], 'LineWidth', 3); 
p2 = polarplot(hdBinsPlot, cconv(hdRatesFast./max(hdRatesFast), win, nHdBins), 'Color', [0.8500 0.3250 0.0980]);
polarplot([hdMeanFast hdMeanFast], [0 hdScoreFast], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 3);

legend([p1 p2], {'Slow', 'Fast'}, 'location', 'best');
title(['Median speed = ' num2str(medSpeed)]);

count = count + 1;
