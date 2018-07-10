subplot(m, n, count);

distVec = linspace(0, boxSize, 10);

meanBorderX = nan(1, 10);
meanBodrerY = nan(1, 10);
% shuffle
for ii_shuffle = 1:1000
    % shuffling and calculations
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps); % create shifted timestamps
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iC = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    % variable of interest
    [~, ~, borderRatesShuffled] = calculate_border_map(iC, vt);
    meanBorderX = mean([meanBorderX; borderRatesShuffled{1}], 'omitnan');
    meanBorderY = mean([meanBodrerY; borderRatesShuffled{2}], 'omitnan');
end

%% x axis
yyaxis right
area(distVec, borderRates{1}, 'FaceAlpha', .5, 'LineStyle', 'none');
ylabel('East-West');
hold on
plot(distVec, meanBorderX);

yyaxis left
area(distVec, borderRates{2}, 'FaceAlpha', .5, 'LineStyle', 'none');
ylabel('North-South');
hold on
plot(distVec, meanBorderY);

xlim([0 96.5]);
xlabel('Distance from Wall (cm)');

legend({'East-West', 'EW Shuffle', 'North-South', 'NS Shuffle'});

count = count + 1;