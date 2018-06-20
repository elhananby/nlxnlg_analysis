dataFolder = fullfile('D:\', 'new_exp_data', 'Cells', '*.mat');
files = subdir(dataFolder);

for iiFile = 1:length(files)
    fileToLoad = files(iiFile).name;
    load(fileToLoad);
    
    vt = cVt;
    clearvars cVt;
    
    vtKeepIdx = index_to_keep(vt, p, s);
    cKeepIdx = index_to_keep(c, p, s);
    
    [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx);
    
    [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);
    
    [speedOccupancy, speedSpikes, speedRates] = calculate_speed_map(c, vt, cKeepIdx, vtKeepIdx);
    
    spikeISI = diff(c.timestamps.*1e-6);

    [rSpikeTrain, lagsSpikeTrain] = xcorr(c.spikeTrain, 500);
end 