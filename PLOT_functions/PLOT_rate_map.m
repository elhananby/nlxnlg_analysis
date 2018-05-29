%% rate map
subplot(m, n, count);
pcolor(posRatesSmooth');
colormap(jet);
shading flat;
axis tight;
axis off;
colorbar;

count = count + 1;
