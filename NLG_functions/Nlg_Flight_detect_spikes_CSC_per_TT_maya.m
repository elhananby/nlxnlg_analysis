function filename_out = Nlg_Flight_detect_spikes_CSC_per_TT_maya(p, ii_tetrode)

dbstop if error

% check if we need to run this function at all
Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');
filename_out = fullfile(p.path_dataout, p.datadir_out,...
    sprintf('TT%i_cleaned', ii_tetrode));

if exist([filename_out '.mat'], 'file')
    fprintf('TT %i already processed - skipping\n', ii_tetrode);
    return;
end

% PARAMETERS:
x_sep_spike_thres = 4; %
amplitude_factor = 1;
r_threshold = p.r_threshold;

library_file = 'library_of_acceptable_spike_shapes.mat'; % File containing the Library Of Acceptable Spike Shapes

Last_Spike_IX = 0; % Initialize
SPK_waveforms = {};% Initialize
SPK_timestamp = [];% Initialize

% Detect crossing-threshold spikes and save their waveforms:
Spike_waveforms_combined = [];
Spike_timestamps_combined = [];

% need to fix it to go over only the channels of the specific tetrode
active_channels = reshape(find(p.active_channels == 1), 4,4)';
active_channels = active_channels(ii_tetrode, :);

fprintf('1)\tSpike detection');
tic
for curr_chan_idx = 1:4
    
    filename_CSC_filtered_in = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered.ncs', active_channels(curr_chan_idx)));
    
    [nlxTimestamps, nlxSamples, Header] = Nlx2MatCSC( filename_CSC_filtered_in, [1 0 0 0 1], 1, 1, []);
    
    Samples{curr_chan_idx} = nlxSamples(:);
    dts = nlxTimestamps(2) - nlxTimestamps(1);
    Timestamps{curr_chan_idx} = linspace(nlxTimestamps(1), nlxTimestamps(end) + dts, length(nlxSamples(:)));
    
    Spike_threshold_uV_units = 4*median(abs(Samples{curr_chan_idx}))./0.6745; % Quiroga 2004
    %     fprintf('Threshold = %.2fuV\n', Spike_threshold_uV_units);
    
    thres_cross_IX{curr_chan_idx} = find(abs(Samples{curr_chan_idx}) > Spike_threshold_uV_units);
    
    % place '1' at threshold crossing:
    thres_cross_vec{curr_chan_idx} = abs(Samples{curr_chan_idx}) > Spike_threshold_uV_units;
    
    % Find the start and end segment of each spiking event on each channel of the tetrode:
    diff_thres_cross_vec{curr_chan_idx} = diff(thres_cross_vec{curr_chan_idx});
    SPK_start{curr_chan_idx} = find(diff_thres_cross_vec{curr_chan_idx} == 1)+1;
    SPK_End{curr_chan_idx} = find(diff_thres_cross_vec{curr_chan_idx} ==(-1) )+1;
    
    %Check for the case that we start with a negetive spike phase
    if (length(SPK_End{curr_chan_idx})>length(SPK_start{curr_chan_idx}))
        SPK_End{curr_chan_idx}(1) = [];
    end
    
    %Check for the case that we start with a negetive spike phase
    if ~isempty(SPK_End{curr_chan_idx})
        if SPK_start{curr_chan_idx}> SPK_End{curr_chan_idx}(1)
            SPK_End{curr_chan_idx}(1) = [];
            SPK_start{curr_chan_idx}(end) = [];
        end
    end
    
    % Extract local segments and find the maxima of each for each channel.
    temp_SPK_events = []; % Initialize
    SPK_max{curr_chan_idx} = [];
    counter = 0;
    SPK_IX{curr_chan_idx} = [];
    
    for jj = 1:length(SPK_End{curr_chan_idx})
        counter = counter + 1;
        temp_SPK_events = Samples{curr_chan_idx}(SPK_start{curr_chan_idx}(jj):SPK_End{curr_chan_idx}(jj));
        [temp_SPK_max,IX] = max(abs(temp_SPK_events));
        temp_IX_vec = SPK_start{curr_chan_idx}(jj):SPK_End{curr_chan_idx}(jj);
        max_IX = temp_IX_vec(IX);
        SPK_IX{curr_chan_idx}(1,jj) = max_IX; % Correct for the real timestamp relative to the entire recording session
    end
    
    %         SPK_IX{curr_chan_idx} = [];
end

% Do coincedence detection to identify single spikes:
% First find spikes which are too close apart to be considered as seperate spikes,
% By two close I mean to have a temporal difference of < x_sep_spike_thres bins between them.
SPK_All_Ch_combined = SPK_IX{1};
fprintf('\t%.2f Sec\n', toc);

fprintf('2)\tCoincidence detection\n');
outTic = tic;

for curr_chan_idx = 2:4
    % First loop over the spikes detected on the first channel and compare their seperation from those detected on the second channel
    % Save only those which are seperated by a minimal number of bins to avoid counting the same spike twice
    fprintf('Channel %i ', curr_chan_idx);
    inTic = tic;
    
    for ii_spike = 1:length(SPK_All_Ch_combined)
        line = fprintf('%.2f%% %i/%i', ii_spike*100/length(SPK_All_Ch_combined), ii_spike, length(SPK_All_Ch_combined));
        shared_IXs = [];
        current_spike_IX = SPK_All_Ch_combined(ii_spike);
        shared_IXs = find(abs(SPK_IX{curr_chan_idx} - current_spike_IX)<= x_sep_spike_thres);
        if ~isempty(shared_IXs) % i.e., the same spike is detected twice
            SPK_IX{curr_chan_idx}(shared_IXs) = [];
        end
        fprintf(repmat('\b', 1, line));
    end
    
    fprintf('\t%.2f Sec\n', toc(inTic));
    % Now that we do not have the same spikes on two channels we can merge the spike IXs of both channels,
    % as those represent different spikes
    SPK_All_Ch_combined = unique([SPK_All_Ch_combined, SPK_IX{curr_chan_idx}]);
end

% sort them to preserve their temporal order:
SPK_IXs_All_Ch_combined_sorted = sort(SPK_All_Ch_combined);
% SPK_timestamp = Timestamps{1}(SPK_IXs_All_Ch_combined_sorted);
fprintf('\t%.2f Sec\n', toc(outTic));

%===============================================================
% Extract the spike waveforms on all 4 channels of the tetrode:
%===============================================================
% Extract the waveform and timestamps of each spike (as in Neuralynx, the peak
% will be the 8th sample out of over all 32 samples of
% the spike shape vec.

fprintf('3)\tExtract spike waveform');
tic

current_file_spike_counter = 0;
for ii_spike = 1:length(SPK_IXs_All_Ch_combined_sorted)
    current_spike_max_IX = SPK_IXs_All_Ch_combined_sorted(ii_spike);
    
    if (current_spike_max_IX+24 <= length(Timestamps{1})...
            && (current_spike_max_IX - 7 > 0)) % i.e., we are NOT cutting the spike in the middle
        current_file_spike_counter = current_file_spike_counter + 1;
        SPK_timestamp(1, Last_Spike_IX + current_file_spike_counter) = Timestamps{1}(current_spike_max_IX);
        
        for curr_chan_idx = 1:4
            if p.active_channels(find(p.TT==ii_tetrode),curr_chan_idx)==1 %i.e. - valid channel
                
                SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter} = Samples{curr_chan_idx}(current_spike_max_IX-7:1:current_spike_max_IX+24);
                %plot(SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter},colors(curr_chan_idx))
                [temp_max(curr_chan_idx,1) temp_max(curr_chan_idx,2)] = max(SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter});
            else
                SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter} = zeros(1,32);
            end
        end
        
    end % if
end % for - ii_spike
Last_Spike_IX = Last_Spike_IX + current_file_spike_counter;
fprintf('\t%.2f Sec\n', toc);


%% Clean artifacts using the library of acceptable spike shapes - All spike waveforms:
fprintf('4)\tClean artifacts\n');
tic

load(library_file); % Load the library of acceptable spike shapes;

vector_of_accepted_spikes = nan( 1, length(SPK_waveforms) ) ; % Initialize
vector_of_max_r_values = nan( 1, 1 ) ;

parfor ii_spike = 1:size(SPK_waveforms, 2)
    spike_shape_4channels = zeros(32,4);
    
    for jj_channel = 1:size(SPK_waveforms, 1)
        spike_shape_4channels(:,jj_channel) = SPK_waveforms{jj_channel, ii_spike};
    end
    
    % Choose the channel # for which the spike has the largest height:
    [ ~,  idx_channel_max_height ] = max( abs(spike_shape_4channels(8,:))) ;%GILY - changed to max of al channels just in the 8th spot (that was defined as peak)

    spike_shape = spike_shape_4channels( :, idx_channel_max_height )' ;

    if (std( spike_shape(2:end-1) ) == 0 ) % If this is a completely FLAT "spike", I cannot compute CORRCOEF, so I will set r = 0
        vector_of_max_r_values( ii_spike ) = 0 ;  % Set r = 0 in this case
    else % If this spike DOES have some shape (this is the case basically for ALL the recorded waveforms)
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag 0 :
        xxx_lags = 2 : 31 ; % Use CENTRAL 30 points
        ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
        rrr_vec_lag_0 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag (+1) :
        xxx_lags = 1 : 30 ; % Use FIRST 30 points (RIGHT shift of the "Acceptable Spike Shapes matrix")
        ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
        rrr_vec_lag_plus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag (-1) :
        xxx_lags = 3 : 32 ; % Use LAST 30 points (LEFT shift of the "Acceptable Spike Shapes matrix")
        ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
        rrr_vec_lag_minus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Save the MAXIMAL r value -- use the maximal correlation among the 3 lags (-1,0,+1):
        vector_of_max_r_values( ii_spike ) = max( [ rrr_vec_lag_0  rrr_vec_lag_plus1  rrr_vec_lag_minus1 ] );
    end

    % Determine if this spike should be Accepted (set value to '1') or Rejected (set value to '0'):
    vector_of_accepted_spikes( ii_spike ) = ( vector_of_max_r_values( ii_spike )  >=  r_threshold );
    
end % End "Loop over spikes extracted from the Ntt file"

% Find the IXs of the accepted spike waveforms and store the accepted and
% not-accepted waveforms speratly.
IX_accepted = find(vector_of_accepted_spikes == 1);
IX_NO_accepted = find(vector_of_accepted_spikes == 0);

% %  if we want to skip comparison to library
%     vector_of_accepted_spikes = ones(1, length(SPK_waveforms));
%     IX_accepted = 1:1:length(SPK_waveforms);
%     IX_NO_accepted = [];

fprintf('\t%.2f Sec\n', toc);

% Extract the accepted (and not accepted) spikes and define new varialbes:
fprintf('5)\tSpike extraction');
tic
Spike_waveforms_accepted = cell(4, numel(IX_accepted));
Spike_waveforms_NO_accepted = cell(4, numel(IX_NO_accepted));

counter_accepted = 0;
counter_NO_accepted = 0;
for ii = 1:length(vector_of_accepted_spikes)
    if vector_of_accepted_spikes(ii) == 1
        counter_accepted = counter_accepted + 1;
        
        for curr_chan_idx=1:4
            Spike_waveforms_accepted{curr_chan_idx,counter_accepted} = SPK_waveforms{curr_chan_idx, :};
        end
    else
        counter_NO_accepted = counter_NO_accepted + 1;
        for curr_chan_idx=1:4
            Spike_waveforms_NO_accepted{curr_chan_idx,counter_NO_accepted} = SPK_waveforms{curr_chan_idx, :};
        end
    end
end
fprintf('\t%.2f Sec\n', toc);

%---------------------------------------------------------------
% Convert into NSE file for spike sorting - All spike waveforms
%---------------------------------------------------------------
fprintf('6)\tConvert and save');
tic
Timestamps_accepted_spikes = SPK_timestamp(IX_accepted);

%change Spike_waveforms_accepted to be a 32 x 4 x num_records matrix.
%[numRecs, numPoints] = size(Spike_waveforms_accepted);
[numCh, numRec] = size(Spike_waveforms_accepted);
spikes = zeros(32,4,numRec);

for rec=1:numRec
    %rec/numRec
    for channel = 1:numCh
        current_channel_waveform = Spike_waveforms_accepted{channel,rec};
        %for point=1:32
        %             spikes(point, channel, rec) = current_channel_waveform(point)*amplitude_factor;
        %end
        spikes(:, channel, rec) = current_channel_waveform*amplitude_factor;
    end
end

Timestamps_accepted_spikes_per_TT = Timestamps_accepted_spikes;
spikes_per_TT = spikes;

save(filename_out, 'spikes_per_TT','Timestamps_accepted_spikes_per_TT')
fprintf('\t%.2f Sec\n', toc);

end
