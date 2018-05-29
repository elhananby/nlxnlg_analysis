function NLG_clean_artifacts(p)
warning('off','all')
%-----------------------------------------------------------------------
% PARAMETERS:
PLT = 0; % '1' if you we want to plot for DBG or '0' otherwise.
AWidth = 33; % Minimal interval (in samples) between two adjacent artifact-threshold-crossings:
% two threshold-crossings with smaller interval will be MERGED
% into ONE artifact. The reason to choose this number
% (equvalent to 2 ms) is that spikes close than 1 ms to
% artifact start or end are removed anyway later on.
min_artifact_length = 4; % Minimum length of artifact (assuming peak of spikes) -- spikes are ALWAYS SHORTER than this.
buff = 33; % buffer around artifacts ~ 1ms
Confidence_boundary_addition = 130; % Artifact threshold - as percentage of the maximal and minimal values in the crudly sorted sleep data
pos_threshold = 1000; % positive threshold for artifact removal (in uv)
neg_threshold = -1000; % positive threshold for artifact removal (in uv)

Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');
flag = zeros(1,16);

%% Find and extract artifacts from the session blocks for each channel
for ii_channel = 1:16 % temporary
    
    % load data and setup
    tic
    fprintf('Cleaning CSC of Day %i from Channel %i', p.day, ii_channel);
    
    filename_out = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered_cleaned', ii_channel));
    
    if exist([filename_out '.ncs'], 'file')
        flag(ii_channel) = 1;
        fprintf('\tAlready processed, skipping %.2f\n', toc);
        continue;
    end
    
    filename_CSC_Filtered_in = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered.ncs', ii_channel));
    
    [Timestamps, Samples] = Nlx2MatCSC( filename_CSC_Filtered_in, [1 0 0 0 1], 0, 1, []);
    
    % artifact cleaning
    spreadSamples = abs(Samples(:));
    
    [~, locs, w, ~] = findpeaks(spreadSamples,...
        'MinPeakHeight', 500,...
        'WidthReference', 'halfprom');
    peakseek
    
    for ii_locs = 1:length(locs)
        idxToClean = locs(ii_locs) - ceil(w(ii_locs)*2) : locs(ii_locs) + ceil(w(ii_locs)*2);
        
        belowZero = idxToClean <= 0;
        idxToClean(belowZero) = 1;
        
        aboveMax = idxToClean > length(spreadSamples);
        idxToClean(aboveMax) = length(spreadSamples);
        
        idxToClean = unique(idxToClean);
        
        totalIdxToClean{ii_channel, ii_locs} = idxToClean;
        
    end % locs
    
    fprintf('\t%.2f Seconds\n', toc);
    
end % channels

for ii_channel = 1:16
    
    if flag(ii_channel) == 1
        continue;
    end
    
    filename_CSC_Filtered_in = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered.ncs', ii_channel));
    
    filename_out = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered_cleaned', ii_channel));
    
    [Timestamps, Samples, Header] = Nlx2MatCSC( filename_CSC_Filtered_in, [1 0 0 0 1], 1, 1, []);
    dt = mean(diff(Timestamps))/512;
    idxToRemove = unique([totalIdxToClean{:}]);
    
    spreadSamples = Samples(:);
    spreadTimestamps = linspace(Timestamps(1),...
        Timestamps(end) + mean(diff(Timestamps)) - dt,... 
        length(spreadSamples));
    
    spreadSamples(idxToRemove) = [];
    spreadTimestamps(idxToRemove) = [];
    
    leftOver = mod(length(spreadTimestamps), 512);
    
    spreadSamples(end-leftOver+1:end) = [];
    spreadTimestamps(end-leftOver+1:end) = [];
    
    reshapeSamples = reshape(spreadSamples,...
        [512, length(spreadSamples)/512]);
    
    reshapeTimestamps = reshape(spreadTimestamps,...
        [512, length(spreadTimestamps)/512]);
    
    Samples = reshapeSamples;
    Timestamps = reshapeTimestamps(1, :);
    
    % save ncs
    AppendToFileFlag = 0;
    ExportMode = 1;
    ExportModeVector = [];
    FieldSelectionFlags = [1 0 0 0 1 1];
    
    Mat2NlxCSC([filename_out '.ncs'], AppendToFileFlag, ExportMode, ExportModeVector,...
        FieldSelectionFlags, Timestamps, Samples, Header);
end
