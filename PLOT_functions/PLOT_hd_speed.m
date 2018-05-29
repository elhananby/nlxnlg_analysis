subplot(m, n, count);

speedIdxPos = vt.speed <= 1;
speedIdxSpk = c.speed <= 1;

hdBins = -pi : 2*pi/nHdBins : pi;

hdRatesSlow = histcounts(c.poshd(speedIdxSpk), hdBins) ./ (histcounts(vt.poshd(speedIdxPos), hdBins) .* (dt * 1e-6));
hdRatesFast = histcounts(c.poshd(~speedIdxSpk), hdBins) ./ (histcounts(vt.poshd(~speedIdxPos), hdBins) .* (dt * 1e-6));

hdBins = linspace(-pi, pi, 60);
polarplot(hdBins, smooth(hdRatesSlow));
hold on
polarplot(hdBins, smooth(hdRatesFast));
legend('<= 1 cm/s', '> 1 cm/s', 'location', 'best');

count = count + 1;