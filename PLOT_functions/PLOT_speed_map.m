%% speed
subplot(m ,n, count);
speedBins = linspace(min(vt.speed), max(vt.speed), nSpeedBins-1);
yyaxis left;
area(speedBins, speedOccupancy, 'FaceAlpha', .5, 'LineStyle', 'none');
ylabel('Time (sec)');

yyaxis right;
area(speedBins, speedRates, 'FaceAlpha', .5, 'LineStyle', 'none');
xlabel('Speed (cm/s)');
ylabel('Firing Rate (Hz)');

title(sprintf('Speed correlation = %.2f', speedScore));

axis tight
box off
count = count + 1;