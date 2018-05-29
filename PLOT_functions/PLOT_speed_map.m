%% speed
subplot(m ,n, count);
speedBins = linspace(min(vt.speed(vtKeepIdx)), max(vt.speed(vtKeepIdx)), nSpeedBins-1);
yyaxis left;
area(speedBins, speedOccupancy, 'FaceAlpha', .5, 'LineStyle', 'none');
ylabel('Time (sec)');

yyaxis right;
area(speedBins, speedRates, 'FaceAlpha', .5, 'LineStyle', 'none');
xlabel('Speed (cm/s)');
ylabel('Firing Rate (Hz)');

axis tight
box off
count = count + 1;