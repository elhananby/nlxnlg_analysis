subplot(m, n, count);

speedIdxPos = vt.speed <= 1;
speedIdxSpk = c.speed <= 1;

hdBins = -pi : 2*pi/nHdBins : pi;

hdRatesSlow = histcounts(c.poshd(speedIdxSpk), hdBins) ./ (histcounts(vt.poshd(speedIdxPos), hdBins) .* dt );
hdRatesFast = histcounts(c.poshd(~speedIdxSpk), hdBins) ./ (histcounts(vt.poshd(~speedIdxPos), hdBins) .* dt);

hdBins = linspace(-pi, pi, 60);
polarplot(hdBins, cconv(hdRatesSlow./max(hdRatesSlow), win, nHdBins));
hold on
polarplot(hdBins, cconv(hdRatesFast./max(hdRatesFast), win, nHdBins));

legend('<= 1 cm/s', '> 1 cm/s', 'location', 'best');

count = count + 1;