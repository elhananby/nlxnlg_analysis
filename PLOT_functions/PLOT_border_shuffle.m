%% shuffle border score
subplot(m, n, count);

shuffledBorderScore = zeros(1, 1000);

for ii_shuffle = 1:1000
    
    % shuffling and calculations
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    % calculate variable
    [~, ~, ~, posRatesSmooth] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx); % calculate rate map
    shuffledBorderScore(ii_shuffle) = calculate_border_score(posRatesSmooth); % calculate border score
end

borderThr = quantile(shuffledBorderScore, 0.99);
histogram(shuffledBorderScore, 40);
hold on;
line([borderThr, borderThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
line([borderScore, borderScore], ylim, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');

title('Shuffled Border Scores');
xlabel('Border Score');
ylabel('Count');

count = count + 1;
