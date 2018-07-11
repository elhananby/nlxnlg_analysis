% head direction shuffling
subplot(m, n, count);

shuffledSpeedScore = zeros(1, 1000);
for ii_shuffle = 1:1000
    
    % shuffling and calculations
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    % variable of interest

    [~, ~, ~, shuffledSpeedScore(ii_shuffle)] = calculate_speed_map(iC, vt);
end

speedThr = quantile(shuffledSpeedScore, 0.99);
histogram(shuffledSpeedScore, 40);
hold on;
line([speedThr, speedThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
line([speedScore speedScore], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');

title('Shuffled Speed Scores');
xlabel('Speed Score');
ylabel('Count');

count = count + 1;