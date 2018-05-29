% head direction shuffling
subplot(m, n, count);

shuffledHdScore = zeros(1, 1000);
for ii_shuffle = 1:1000
    
    % shuffling and calculations
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    % variable of interest
    [~, ~, shuffledHdRates, ~] = calculate_hd_map(iC, vt, iCKeepIdx, vtKeepIdx);
    hdBins = linspace(-pi, pi, nHdBins);
    shuffledHdScore(ii_shuffle) = circ_r(hdBins', shuffledHdRates');
end

hdThr = quantile(shuffledHdScore, 0.99);
histogram(shuffledHdScore, 40);
hold on;
line([hdThr, hdThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
hdScore = circ_r(hdBins', hdRates');
line([hdScore hdScore], ylim, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');

title('Shuffled HD Scores');
xlabel('HD Score');
ylabel('Count');

count = count + 1;