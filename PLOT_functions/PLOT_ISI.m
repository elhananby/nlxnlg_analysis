%% ISI
subplot(m, n, count);
x = logspace(-4, 1, 100);
histogram(spikeISI, x);
set(gca, 'xscale', 'log');
hold on
line([0.002 0.002], [0 max(histcounts(log(spikeISI), 100))], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('ISI Count');
title(sprintf('%i# Spikes\n%.2fHz',...
    length(c.timestamps),...
    length(c.timestamps)/((vt.timestamps(end)-vt.timestamps(1))*1e-6)));
box off
count = count + 1;