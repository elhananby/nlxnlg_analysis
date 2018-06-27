% loop over variables of interest
% randVec = rand(1, 1000);
for jj = 1:1000
    
    % shuffling and calculations
%     shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps, randVec(jj)); % create shifted timestamps
    shiftedTimestamps = permute_timestamps(c.timestamps);
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    %% head direction
    % total
    [~, ~, shuffleHdRates, ~] = calculate_hd_map(iC, vt, iCKeepIdx, vtKeepIdx);
    shuffledHdScore(jj) = circ_r(hdBins', shuffleHdRates');
    
    % first half
    halfToKeep{1} = find(vt.timestamps <= vt.timestamps(ceil(length(vt.timestamps)/2)))';
    halfToKeep{2} = find(c.timestamps <= vt.timestamps(ceil(length(vt.timestamps)/2)))';
    
    [~, ~, shuffleHdRates1, ~] = calculate_hd_map(iC, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    shuffledHdScore1(jj) = circ_r(hdBins', shuffleHdRates1');
    
    % second half
    halfToKeep{1} = find(vt.timestamps > vt.timestamps(ceil(length(vt.timestamps)/2)))';
    halfToKeep{2} = find(c.timestamps > vt.timestamps(ceil(length(vt.timestamps)/2)))';
    
    [~, ~, shuffleHdRates2, ~] = calculate_hd_map(iC, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    shuffledHdScore2(jj) = circ_r(hdBins', shuffleHdRates2');
    
    % slow
    halfToKeep{1} = find(vt.speed <= 1)';
    halfToKeep{2} = find(c.speed <= 1)';
    
    [~, ~, shuffleHdRatesSlow, ~] = calculate_hd_map(iC, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    shuffledHdScoreSlow(jj) = circ_r(hdBins', shuffleHdRatesSlow');
    
    % fast
    halfToKeep{1} = find(vt.speed > 1)';
    halfToKeep{2} = find(c.speed > 1)';
    
    [~, ~, shuffleHdRatesFast, ~] = calculate_hd_map(iC, vt,...
        intersect(cKeepIdx, halfToKeep{2}),...
        intersect(vtKeepIdx, halfToKeep{1}));
    
    shuffledHdScoreFast(jj) = circ_r(hdBins', shuffleHdRatesFast');

    %     % border clean
%     [~, ~, ~, posRatesSmooth] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx);
%     try % some border scores cant be calcuated due to shuffling that leaves all zeros - quick fix
%         [shuffledBorderScore(jj), shuffledBorderDir(jj)] = calculate_border_score(posRatesSmooth);
%     catch
%         continue;
%     end
%     
%     % border original
%     [~, ~, ~, ~, posRatesSmoothUnclean] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx);
%     try
%         [shuffledBorderScoreUnclean(jj), shuffledBorderDirUnclean(jj)] = calculate_border_score(posRatesSmoothUnclean);
%     catch
%         continue;
%     end
    
end