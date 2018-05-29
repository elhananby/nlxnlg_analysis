%% spike shape
subplot(m, n, count);

% find channel with maximum value
[~, maxChan] = max(max(mean(c.spikeShape, 3)));

% plot all spike lines
plot(squeeze(c.spikeShape(:, maxChan, :)), 'Color', [220 220 220 50]./255, 'LineWidth', .5);
hold on;

% plot only mean of spike
plot(mean(squeeze(c.spikeShape(:, maxChan, :)), 2), 'k', 'LineWidth', 1);

xlabel('Time (ms)');
ylabel('Voltage (uV)');
box off;

count = count + 1;