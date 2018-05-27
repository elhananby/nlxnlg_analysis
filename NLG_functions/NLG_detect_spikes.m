function NLG_detect_spikes(p)
<<<<<<< HEAD
dbstop if error
filename_out_temp = fullfile(p.path_dataout, p.datadir_out, 'TT1_cleaned.ntt');
if exist(filename_out_temp, 'file')
    return;
end

fprintf('Spike Extraction:\n');

tts = reshape([1:16],4,4);

for ii_tetrode = 1:length(param.tetrodes.use_for_sorting)

    filename = cell(1,4);
    Spiketimeidx = cell(1,4);
    csc = cell(4,1);
    path = [param.path.CSC_spikes,'spikes_bat',param.bat,'_',param.day,'_TT'];
    tetrodes_use_for_sorting = param.tetrodes.use_for_sorting;
    x_sep_spike_thres = param.spikes.x_sep_spike_thres;
    thresh = zeros(1,4);
   
    % detect spikes
    parfor ch = 1:4
        
        filename = [path,num2str(tts(ch,tetrodes_use_for_sorting(ii_tetrode)))];
        csc{ch,1} = load(filename);
        thresh(ch) = 2.8 * median(abs(csc{ch,1}.Spikes.data.Samples(:))/0.6745);
        [spikeidx, spikepk]  = peakseek(csc{ch,1}.Spikes.data.Samples, x_sep_spike_thres, thresh(ch));
        [spikeidx, ii] = sort(spikeidx);
        spikepk = spikepk(ii);
        spikeidx = spikeidx((spikeidx>7 & spikeidx<=length(csc{ch,1}.Spikes.data.Samples)-24));
        spiketime = csc{ch,1}.Spikes.data.Timestamps(spikeidx);
        spiketimeidx{ch,1} = [spiketime(:)';spikeidx(:)'];  
        
    end
    
    threshold(:,ii_tetrode) = thresh(:); 
    SamplingFreq = csc{ii_tetrode}.Spikes.params.CSC_SamplingFreq;
    
    fprintf('Processing Tetrode %i...\n', ii_tetrode);
    
    %% coincedent detection (matrixwise implementation) and removal of spikes which fall within a refractory period
    startidx = min([spiketimeidx{1}(2,:)' ;spiketimeidx{2}(2,:)';spiketimeidx{3}(2,:)';spiketimeidx{4}(2,:)']);
    endidx  =  max([spiketimeidx{1}(2,:)' ;spiketimeidx{2}(2,:)';spiketimeidx{3}(2,:)';spiketimeidx{4}(2,:)']);
    len = startidx - endidx +1;
    M = zeros(4,len);
    
    for i = 1:4
        M(i,spiketimeidx{i}(2,:)) = csc{i,1}.Spikes.data.Samples(spiketimeidx{i}(2,:));
    end
    
    M = max(M);
    
    spikeidx = peakseek(M,22,min(thresh(:)));
    spiketime = csc{1,1}.Spikes.data.Timestamps(spikeidx);
    clear M;
    clear xx;
    clear idx;

    %% calculate spikewave
    spikewave = zeros(32,size(spiketime,2),4);
    mat= repmat([-7:24]',1,length(spikeidx));
    spike_idx = repmat(spikeidx,[32,1]);
    spike_idx = spike_idx+mat;
    parfor ch=1:4
        spikewave(:,:,ch) = csc{ch,1}.Spikes.data.Samples(spike_idx);
    end
    clear csc;
    
    %% compare to acceptable spike shapes
    eval(['load ', param.path.spike_shapes_lib]); % Load the library of acceptable spike shapes;
    [~,iimax ] = max(squeeze(max(spikewave)),[],2);
    waves = zeros(size(spikewave,1),size(spikewave,2));
    for i=1:size(waves,2)
        waves(:,i) = spikewave(:,i,iimax(i));
    end
    library_of_acceptable_spike_shapes = library_of_acceptable_spike_shapes';
    c = zeros(size(library_of_acceptable_spike_shapes,2),size(waves,2));
    parfor ii=1:size(library_of_acceptable_spike_shapes,2)
        mat = repmat(library_of_acceptable_spike_shapes(:,ii),1,size(waves,2));
        c1 = corrcoeff(waves(2:end-1,:),mat(2:end-1,:));
        c2 = corrcoeff(waves(1:end-2,:),mat(2:end-1,:));
        c3 = corrcoeff(waves(3:end,:),mat(2:end-1,:));
        c(ii,:) = max([c1;c2;c3]);
    end
    c = max(c);
    acceptedidx = find(c>=0.8);
    Timestamps_accepted_spikes_TT{ii_tetrode}=spiketime(acceptedidx);
    spikes_TT{ii_tetrode}=spikewave(:,acceptedidx,:);
    clear c;
end

%%  Coincidence-Detection across Tetrodes to eliminate artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(param.tetrodes.use_for_sorting)>1
    idx_coincidence_vec = cell(length(param.tetrodes.use_for_sorting),1);
    
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}) % Loop over the spikes of the FIRST tetrode
        temp_stack = cell(length(param.tetrodes.use_for_sorting),1);
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        temp_stack{1,1} = ii_spikes;
        for TT=2:length(param.tetrodes.use_for_sorting)
            temp_stack{TT,1} = find( abs( Timestamps_accepted_spikes_TT{TT} - t_spike ) <= param.spikes.coincidence_window ); % THE COINCIDENCE DETECTION
        end
        ii = ~cellfun(@isempty,temp_stack);
        if sum(ii)>=length(param.tetrodes.use_for_sorting)
            for jj=1:length(param.tetrodes.use_for_sorting)
                idx_coincidence_vec{jj,1} = cat(2,idx_coincidence_vec{jj,1}, temp_stack{jj,1});
            end
        end
    end
end
   

for i=1:length(param.tetrodes.use_for_sorting)
    idx_coincidence_vec{i} = unique(idx_coincidence_vec{i});
    if ~isempty(idx_coincidence_vec{i})
        Timestamps_accepted_spikes_TT{i}(idx_coincidence_vec{i})=[];
        spikes_TT{i}(:,idx_coincidence_vec{i},:)=[];
    end
end


%% Save the data in Neuralynx NTT files for three cases:
%%%%%%%%%%%%%%%%%%%%%%%%%

for ii_Tetrode=1: length(param.tetrodes.use_for_sorting)   
          
    base_output_dir						= [param.path.analysis, 'preprocessing_spikes4',filesep];
    if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end
    filename_out = [base_output_dir, 'spikes_',param.DirRaw,'_TT', num2str(param.tetrodes.use_for_sorting(ii_Tetrode))];
    
    % constract header for NTT file
    load(param.path.NTT_header);
    l = ['-ThreshVal ',num2str(threshold(1,ii_Tetrode)),' ',...
                    num2str(threshold(2,ii_Tetrode)),' ',...
                    num2str(threshold(3,ii_Tetrode)),' ',...
                    num2str(threshold(4,ii_Tetrode))]; 
    Header{end+1} = l;
    Header = MakeNlxHeader(Header);
    %Header = Header'; 
    
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{ii_Tetrode};
    spikes=spikes_TT{ii_Tetrode};
    spikes = permute(spikes,[1,3,2]);
    
    %only export timestamps and data points - Full session
    FieldSelection = [1 0 0 0 1 1];
    
    if ispc
        Mat2NlxSpike([filename_out,'.NTT'],0,1,[],FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    elseif isunix
        if exist([filename_out,'.NTT'],'file'), delete([filename_out,'.NTT']),end
        Mat2NlxTT([filename_out,'.NTT'],0,1,1,length(Timestamps_accepted_spikes(:)),FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    end
end
end
=======

channelMatrix = p.active_channels .* [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
channelIdx = 1;

for iiTetrode = p.TT
    fprintf('\nProcessing Tetrode %i', iiTetrode);
    
    currChannels = nonzeros(channelMatrix(iiTetrode, :))';
    
    %% 1) spike detection using peakseek
    fprintf('\n1) Finding peaks...');
    peakFind = tic;
    
    for iiChannel = 1:4
        
        % get input filename
        filename_CSC_filtered_in = fullfile(p.path_dataout, p.datadir_out, 'nlx_data',...
            sprintf('CSC%i_filtered.ncs', channelIdx));

        % read using nlx
        [nlxTimestamps, nlxSamples, Header] = Nlx2MatCSC( filename_CSC_filtered_in, [1 0 0 0 1], 1, 1, []);
        
        % this calculates the maximum seperation between samples based on
        % finding the sampling rate and calculating the length of one spike
        % (1000 ms)
%         seperationIdx = ceil(1000/ ... % spike width in microseconds (1000microseconds = 1 milliseconds)
%             (1/str2double(cell2mat(regexpi(Header{8}, '(\d*)\.(\d*)', 'match'))) * 1e6))/4; % sampling rate in microseconds

        seperationIdx = 4;

        Samples{iiChannel} = nlxSamples(:); % get spread samples
        
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
    for iiChannel = 1:4
        
        % in each spike index for each channel, put the measure sample
        % value from the raw data
        M(iiChannel, spikesIdx{iiChannel}) = Samples{iiChannel}(spikesIdx{iiChannel});
    end
    
    % get the maximum value from each column (meaning, get the channel with
    % maximum value)
    M = max(M);
    
    % use peakseek again to find events with a specific seperation between
    % them (22 points?)
    spikesIdx = peakseek(M, 22, min(cell2mat(spikeThreshold)));
    
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
    
    for iiChannel = 1:4
        % use the index matrix to get spike shape
        spikesWaveforms(:, :, iiChannel) = Samples{iiChannel}(spikesWaveformsIdx); 
    end
    
    fprintf('\t%.2f Seconds\n', toc(waveformExtraction));
    
    %% 4) Correlation with library of existing spikes
    fprintf('4) Cleaning artifacts... ');
    artifactCleaning = tic;
    load('D:\Scripts\nlx_analysis\library_of_acceptable_spike_shapes.mat'); % load library file
    
    % find the channel with the maximum spike height for all waveforms
    [ ~, maxChannel ] = max(squeeze(max(spikesWaveforms)), [], 2);
    
    % initialize spike waveforms to test array
    spikeWavesToTest = zeros(size(spikesWaveforms, 1), size(spikesWaveforms, 2));
    
    % and fill it with all spike waveforms for each maximum channel
    for ii = 1:size(spikeWavesToTest, 2)
        spikeWavesToTest(:, ii) = spikesWaveforms(:, ii, maxChannel(ii));
    end
    
    lib = library_of_acceptable_spike_shapes'; % transpose for later
    indexToKeep = []; % initalize
    
    % windowing procedure so we could use pdist2 and not run out of memory
    windowStart = 1;
    windowEnd = 100;
    windowSize = 100;
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
    CellNumbers = zeros(1, length(timestampsToSave{iiTetrode}));
     
    Mat2NlxSpike( sprintf('TT%i.ntt', iiTetrode), 0, 1, [],...
        [1 0 1 0 1 0],  timestampsToSave{iiTetrode}, CellNumbers,...
        spikesToSave{iiTetrode});
    
    clearvars -except timestampsToSave spikesToSave p iiTetrode channelMatrix channelIdx
end

endstop = 1;

end
>>>>>>> Finished working on new spike detection algorithm, still need to debug it and test it against origianl scripts.
