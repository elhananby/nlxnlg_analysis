%% head direction time plot
subplot(m, n, count);

angleDiff = diff(vt.poshd);
correctedPosHd = vt.poshd;
correctedPosHd(abs(angleDiff) >= 2*pi/1.8) = NaN;
plot(vt.timestamps, correctedPosHd, '-k', c.timestamps, c.poshd, '.r');
axis tight
xlim([min(vt.timestamps) max(vt.timestamps)]);
ylim([-pi pi]);

count = count + 1;