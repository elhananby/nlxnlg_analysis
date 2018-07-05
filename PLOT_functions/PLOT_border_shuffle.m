%% shuffle border score
subplot(m, n, count);

shuffledCurves = zeros(size(borderCurve, 2), 1000);

for iiShuff = 1:1000
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps);
    iC = interpolate_shifted_values(shiftedTimestamps, vt);
    
    [~, curve] = calculate_border_score(iC);
    % we are using the original borderIdx variable because we want the
    % shuffling of the original wall
    shuffledCurves(:, iiShuff) = curve(borderIdx, :);
end
meanShuffleCurve = mean(shuffledCurves, 2);
semShuffleCurve = std(shuffledCurves, 0, 2);

distVec = 1:1:96;
boundedline(distVec, borderCurve(borderIdx, :), 0, '-b', distVec, meanShuffleCurve, semShuffleCurve, '--r');
legend('Original', 'Shuffled');
xlabel('Distance from wall (cm)')
ylabel('Proportion of spikes');

count = count + 1;

%% old version
% shuffledBorderScore = zeros(1, 1000);
% 
% for ii_shuffle = 1:1000
%     
%     % shuffling and calculations
%     shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
%     iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
%     iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
%     
%     % calculate variable
%     [~, ~, ~, posRatesSmooth] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx); % calculate rate map
%     shuffledBorderScore(ii_shuffle) = calculate_border_score(posRatesSmooth); % calculate border score
% end
% 
% borderThr = quantile(shuffledBorderScore, 0.99);
% histogram(shuffledBorderScore, 40);
% hold on;
% line([borderThr, borderThr], ylim, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
% line([borderScore, borderScore], ylim, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
% 
% title('Shuffled Border Scores');
% xlabel('Border Score');
% ylabel('Count');
