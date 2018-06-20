hdcells = [];
nHdBins = 60;
line = 0;
tic
for ii = 1:length(cells)
    fprintf(repmat('\b', 1, line));
    line = fprintf('%i/%i %.2fSec %.2f%%', ii, length(cells), toc, ii*100/length(cells));
    
    load(cells{ii});
    % correct struct names for analysis
    c = c;
    vt = cVt;
    
    p = p;
    s = p.S(c.session);
    
    vtKeepIdx = index_to_keep(vt, p, s);
    cKeepIdx = index_to_keep(c, p, s);
    [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);
    
    shuffledHdScore = zeros(1, 1000);
    parfor ii_shuffle = 1:1000
        
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
    if hdScore >= hdThr &&...
            ((s.end_time - s.start_time).*1e-6)/60 >= 10 &&...
            length(c.timestamps)>=300
        hdcells(end+1) = c.cell_number;
    end
end
