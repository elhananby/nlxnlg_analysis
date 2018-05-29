%% position map
subplot(m, n, count);
plot(vt.posx_c, vt.posy_c, '-k', c.posx_c, c.posy_c, '.r');
axis tight
axis off

count = count +1;