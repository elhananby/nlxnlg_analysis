% loop over variables of interest
randVec = rand(1, 1000);
for jj = 1:1000
    
    % shuffling and calculations
    shiftedTimestamps = shift_timestamps(c.timestamps, vt.timestamps, randVec(jj)); % create shifted timestamps
    iC = interpolate_shifted_values(shiftedTimestamps, vt); % calculate interpolated shifted values
    iCKeepIdx = index_to_keep(iC, p, s); % calculate indices to keep for shifted values
    
    % border clean
    [~, ~, ~, posRatesSmooth] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx);
    try % some border scores cant be calcuated due to shuffling that leaves all zeros - quick fix
        [shuffledBorderScore(jj), shuffledBorderDir(jj)] = calculate_border_score(posRatesSmooth);
    catch
        continue;
    end
    
    % border original
    [~, ~, ~, ~, posRatesSmoothUnclean] = calculate_rate_map(iC, vt, iCKeepIdx, vtKeepIdx);
    try
        [shuffledBorderScoreUnclean(jj), shuffledBorderDirUnclean(jj)] = calculate_border_score(posRatesSmoothUnclean);
    catch
        continue;
    end
    
    % head direction
    [~, ~, shuffleHdRates, ~] = calculate_hd_map(iC, vt, iCKeepIdx, vtKeepIdx);
    shuffledHdScore(jj) = circ_r(hdBins', shuffleHdRates');
end