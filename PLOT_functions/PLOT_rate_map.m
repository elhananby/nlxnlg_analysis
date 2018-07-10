%% rate map
subplot(m, n, count);

sc(posRatesSmooth, jet);
colorbar;
set(gca, 'YDir', 'normal'); % imagesc/sc flips the x axis for whichever reason

count = count + 1;
