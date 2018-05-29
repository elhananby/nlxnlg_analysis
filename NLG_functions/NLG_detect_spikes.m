function NLG_detect_spikes(p)
global useGPU;
channelIdx = 1;

if useGPU == 1
    gpuInfo = gpuDevice;
    maxArraySize = gpuInfo.AvailableMemory/8;
else
    memInfo = memory;
    maxArraySize = memInfo.MaxPossibleArrayBytes/8;
end

for iiTetrode = p.TT
    outFile = fullfile(p.path_dataout, p.datadir_out, sprintf('TT%i.mat', iiTetrode));
    if exist(outFile, 'file')
        continue;
    end
    
    fprintf('\nProcessing Tetrode %i', iiTetrode);
    
    currChannels = find(p.active_channels(iiTetrode, :));
    
    %% 1) spike detection using peakseek
    fprintf('\n1) Finding peaks...');
    peakFind = tic;
    
    for iiChannel = currChannels
        
        % get input filename
        filename_CSC_filtered_in = fullfile(p.path_dataout, p.datadir_out, 'nlx_data',...
            sprintf('CSC%i_filtered_cleaned.ncs', channelIdx));

        % read using nlx
        [nlxTimestamps, nlxSamples, Header] = Nlx2MatCSC( filename_CSC_filtered_in, [1 0 0 0 1], 1, 1, []);
        
        % this calculates the maximum seperation between samples based on
        % finding the sampling rate and calculating the length of one spike
        % (1000 ms)
%         seperationIdx = ceil(1000/ ... % spike width in microseconds (1000microseconds = 1 milliseconds)
%             (1/str2double(cell2mat(regexpi(Header{8}, '(\d*)\.(\d*)', 'match'))) * 1e6))/4; % sampling rate in microseconds

        seperationIdx = 4;
        if useGPU
            Samples{iiChannel} = gpuArray(nlxSamples(:)); % get spread samples
        else
            Samples{iiChannel} = nlxSamples(:); % get spread samples
        end
        
        % calculate spread timestamps
        dts = nlxTimestamps(2) - nlxTimestamps(1); 
        Timestamps{iiTetrode} = linspace(nlxTimestamps(1), nlxTimestamps(end) + dts, length(nlxSamples(:)));
        
        % calculate spike threshold for each channel
        % (need to window it for more accurate calculation)
        spikeThreshold{iiChannel} = 2.8 * median( abs(Samples{iiChannel}) ./ 0.6745 );
        
        % use peakseek to find peaks (spikes) based on threshold and
        % seperationidx
        spikesIdx{iiChannel} = peakseek(Samples{iiChannel}, seperationIdx, spikeThreshold{iiChannel});
        
        % keep only values that we can actually save their waveforms
        % (meaning with 7 points before and 24 points after)
        idxRange = spikesIdx{iiChannel} > 7 & spikesIdx{iiChannel} <= length(Samples{iiChannel}-24);
        spikesIdx{iiChannel} = spikesIdx{iiChannel}(idxRange);
        
        % get the spike times
        spikesTimestamps{iiChannel} = Timestamps{iiTetrode}(spikesIdx{iiChannel});
        
        channelIdx = channelIdx + 1;
    end
    
    fprintf('\t%.2f Seconds\n', toc(peakFind));
    
    %% 2) Conicidence detection
    fprintf('2) Within-Tetrode coincidence detection...');
    channelCoincidence = tic;
    
    % get min and maximum points of all spiking events
    startIdx = min([spikesIdx{1}'; spikesIdx{2}'; spikesIdx{3}'; spikesIdx{4}']);
    endIdx = max([spikesIdx{1}'; spikesIdx{2}'; spikesIdx{3}'; spikesIdx{4}']);
    
    % create array for each channel
    len = endIdx - startIdx + 1;

    M = zeros(4, len);

    % go over all channels
    for iiChannel = currChannels
        if useGPU
            Samples{iiChannel} = gather(Samples{iiChannel});
        end
        % in each spike index for each channel, put the measure sample
        % value from the raw data
        M(iiChannel, spikesIdx{iiChannel}) = Samples{iiChannel}(spikesIdx{iiChannel});
    end
    
    % get the maximum value from each column (meaning, get the channel with
    % maximum value)
    M = max(M);
    
    % use peakseek again to find events with a specific seperation between
    % them (22 points?)
    if useGPU
        spikeThreshold = cellfun(@gather, spikeThreshold);
        M = gpuArray(M);
    else
        spikeThreshold = cell2mat(spikeThreshold);
    end
    spikesIdx = peakseek(M, 22, min(spikeThreshold));
    
    % merge all indices found across channels
    % and find their timestamps
    allSpikes = unique(spikesIdx);
    allTimestamps = Timestamps{iiTetrode}(spikesIdx);
    
    fprintf('\t%.2f Seconds\n', toc(channelCoincidence));
    
    %% 3) Extract spike waveform
    fprintf('3) Extracting waveforms...');
    waveformExtraction = tic;
    
    % create matrix of 32 (waveform size) by number of spikes detected
    waveformMat = repmat([-7:24]', [1 length(allSpikes)]);
    
    % add the spike indices to the waveform mat, so we get [-7:24] around spike index
    spikesWaveformsIdx = allSpikes + waveformMat; 
    
    % initialize spikeWaveforms array
    spikesWaveforms = zeros(32, size(allSpikes, 2), length(currChannels));
    
    for iiChannel = currChannels
        % use the index matrix to get spike shape
        spikesWaveforms(:, :, iiChannel) = Samples{iiChannel}(spikesWaveformsIdx); 
    end
    
    fprintf('\t%.2f Seconds\n', toc(waveformExtraction));
    
    %% 4) Correlation with library of existing spikes
    fprintf('4) Cleaning artifacts... ');
    artifactCleaning = tic;
    load('D:\Scripts\nlxnlg_analysis\library_of_acceptable_spike_shapes.mat'); % load library file
    
    % find the channel with the maximum spike height for all waveforms
    [ ~, maxChannel ] = max(squeeze(max(spikesWaveforms)), [], 2);
    
    % initialize spike waveforms to test array
    if useGPU == 1
        spikeWavesToTest = zeros(size(spikesWaveforms, 1), size(spikesWaveforms, 2), 'gpuArray');
        lib = gpuArray(library_of_acceptable_spike_shapes'); % transpose for later
    else
        spikeWavesToTest = zeros(size(spikesWaveforms, 1), size(spikesWaveforms, 2));
        lib = library_of_acceptable_spike_shapes'; % transpose for later
    end
    
    % and fill it with all spike waveforms for each maximum channel
    for ii = 1:size(spikeWavesToTest, 2)
        spikeWavesToTest(:, ii) = spikesWaveforms(:, ii, maxChannel(ii));
    end

    indexToKeep = []; % initalize
    
    % calculate maximum window size we can fit in memory
    if useGPU == 1
        gpuInfo = gpuDevice;
        maxArraySize = gpuInfo.AvailableMemory/8;
    else
        memInfo = memory;
        maxArraySize = memInfo.MaxPossibleArrayBytes/8;
    end
    
    corrArraySize = size(spikeWavesToTest, 2)*size(library_of_acceptable_spike_shapes, 1); % this will be the size of matrix we get out of pdist2
    corrArrayInMem = corrArraySize/maxArraySize; % we divide the number of elements in the matrix in the maximum number of elements we can fit in memory
    
    if corrArrayInMem < 1 % if we can fit in memory less than the whole array
        windowSize = floor(size(spikeWavesToTest, 2) * corrArrayInMem); % set window 
    else
        windowSize = size(spikeWavesToTest, 2) ; % the window will be the entire array
    end
      
    % windowing procedure so we could use pdist2 and not run out of memory
    windowStart = 1;
    windowEnd = windowSize;
    endFlag = 0;
    
    % loop over 
    while endFlag ~= 2
        line = fprintf('%.2f%%', windowEnd*100/size(spikeWavesToTest, 2));
        calcIdx = windowStart:windowEnd;
        
        % calculate correlation using 1-pdist2, which gives us correlation
        % matrix of all columns between the two arrays
        d1 = 1-pdist2(spikeWavesToTest(2:end-1, calcIdx)', lib(2:end-1, :)', 'correlation');
        d2 = 1-pdist2(spikeWavesToTest(1:end-2, calcIdx)', lib(2:end-1, :)', 'correlation');
        d3 = 1-pdist2(spikeWavesToTest(3:end, calcIdx)', lib(2:end-1, :)', 'correlation');
        
        % now find unique indices to save by indexing only the rows
        [d1Index, ~] = find(d1 >= p.r_threshold); d1Index = unique(d1Index);
        [d2Index, ~] = find(d2 >= p.r_threshold); d2Index = unique(d2Index);
        [d3Index, ~] = find(d3 >= p.r_threshold); d3Index = unique(d3Index);
        
        % if the array is not empty
        if ~isempty(unique([d1Index; d2Index; d3Index]))
            % correct for windowing indexing and add it to the current array of indices to keep 
            currentIndexToKeep = unique([d1Index; d2Index; d3Index]) + windowStart - 1;
            indexToKeep = unique([indexToKeep; currentIndexToKeep]);
        end
        
        % update window
        windowStart = windowStart + windowSize;
        windowEnd = windowEnd + windowSize;
        if windowEnd > size(spikeWavesToTest, 2) && endFlag == 0
            windowEnd = size(spikeWavesToTest, 2);
            endFlag = 1;
        elseif windowEnd > size(spikeWavesToTest, 2) && endFlag == 1
            endFlag = 2;
        end
        fprintf(repmat('\b', 1, line));
        
    end
    fprintf('\t%.2f Seconds\n', toc(artifactCleaning));
    
    %% save data
    timestampsToSave{iiTetrode} = allTimestamps(indexToKeep);
    spikesToSave{iiTetrode} = permute(spikesWaveforms(:, indexToKeep, :), [1 3 2]);
    
    Timestamps = timestampsToSave{iiTetrode};
    Samples = spikesToSave{iiTetrode};
    save(outFile, 'Timestamps', 'Samples');
    
    gpuDevice(1); % flush GPU memroy
    clearvars -except timestampsToSave spikesToSave p iiTetrode channelMatrix channelIdx useGPU
end

%% between-tetrode coincidece detection
fprintf('Between-tetrode coincidence detection: ');
TTcoincidence = tic;

coincidenceIdx = cell(length(p.TT), 1);

for iiSpikes = 1:length(timestampsToSave{1})
    tempStack = cell(length(p.TT), 1);
    spikeTime = timestampsToSave{1}(iiSpikes);
    tempStack{1, 1} = iiSpikes;
    for iiTetrode = 2:length(p.TT)
        tempStack{iiTetrode, 1} = find( abs(timestampsToSave{iiTetrode} - spikeTime) <= 500);
    end
    ii = ~cellfun(@isempty, tempStack);
    if sum(ii) >= length(p.TT)
        for jj = 1:length(p.TT)
            coincidenceIdx{jj, 1} = cat(2, coincidenceIdx{jj, 1}, tempStack{jj,1});
        end
    end
end

fprintf('\t%.2f Seconds\n', toc(TTcoincidence));

for iiTetrode = p.TT
    coincidenceIdx{iiTetrode} = unique(coincidenceIdx{iiTetrode});
    if ~isempty(coincidenceIdx{iiTetrode})
        timestampsToSave{iiTetrode}(coincidenceIdx{iiTetrode}) = [];
        spikesToSave{iiTetrode}(:, :, coincidenceIdx{iiTetrode}) = [];
    end
    
    % save data
    
    outFile = fullfile(p.path_dataout, p.datadir_out, sprintf('TT%i.ntt', iiTetrode));
    if exist(outFile, 'file')
        continue;
    end
    
    AppendToFileFlag = [];
    ExportMode = 1;
    ExportModeVector = [];
    FieldSelectionFlags = [1 0 1 0 1 0];
    CellNumbers = zeros(1, length(timestampsToSave{iiTetrode}));
    
    Mat2NlxSpike( outFile, AppendToFileFlag, ExportMode, ExportModeVector,...
        FieldSelectionFlags, timestampsToSave{iiTetrode}, CellNumbers, spikesToSave{iiTetrode});
end

end
