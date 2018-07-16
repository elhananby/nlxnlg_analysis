function Nlg_Flight_detect_spikes_CSC_per_TT(p, ii_Tetrode)

% PARAMETERS:
PLT = 1; % '1' if you we want to plot for DBG or '0' otherwise.
t_extract_before_spikes = 0.25; % in ms
t_extract_after_spikes = 0.75; % in ms
x_sep_spike_thres = 4; %
amplitude_factor = 1;
r_threshold = p.r_threshold;
coincidence_window = 500 ; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes

%-----------------------------------------------------------------------

Day					= p.day;
Bat					= p.animal;

used_channels=1:4;

library_file = 'library_of_acceptable_spike_shapes.mat'; % File containing the Library Of Acceptable Spike Shapes

dir_data_in							= [p.path_dataout p.path_year_bat num2str(Day) '\'];
% %  filename_associated_VT_file			= [ p.path_dataout p.path_year_bat num2str(Day) '\' 'VT_flight_extracted_bat' num2str(Bat) '.mat']


%load the number of CSC chunks to analyze
load ([p.path_datain p.path_year_bat 'CSC_extracted\' ,num2str(Day), '\Num_of_analyzed_files.mat']);
session_blocks(1).file_num = Num_of_pre_sleep_analyzed_files;
session_blocks(2).file_num = Num_of_behav_analyzed_files;
session_blocks(3).file_num = Num_of_post_sleep_analyzed_files;
session_blocks(1).session_name = 'pre-sleep';
session_blocks(2).session_name = 'behavioral';
session_blocks(3).session_name = 'post-sleep';

session_blocks(1).file_name = 'pre_sleep';
session_blocks(2).file_name = 'behav';
session_blocks(3).file_name = 'post_sleep';


session_blocks(1).plot_length = 10000;
session_blocks(2).plot_length = 50000;
session_blocks(3).plot_length = 10000;


%--------------------------------------------
% changed for every recording day.  for this batch file only

Spike_threshold_uV_units = p.Spike_threshold_uV_units;

%-----------------------------------------------------------------------

disp('======================');
disp([' Detecting Spikes from CSC. #' num2str(Day) ' from Tetrode ' num2str(ii_Tetrode)]);
disp('======================');

filename_data_in= [p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\CSC_extracted_bat' num2str(Bat) '_Day' num2str(Day) '_TT', num2str(ii_Tetrode)]; % This is the prefix
filename_associated_artifact_file	= [dir_data_in, 'Clean_artifacts\CSC_artifacts_bat',num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode)]; %this is a prefix

Last_Spike_IX = 0; % Initialize
SPK_waveforms = {};% Initialize
SPK_timestamp = [];% Initialize
for curr_session_block = session_blocks
        for ii_file = 1:curr_session_block.file_num
    disp(['Processing min #',num2str(ii_file),' out of ',num2str(curr_session_block.file_num),' ', curr_session_block.session_name, ' files...'])
            
            %Load the current file data
            current_filename = [filename_data_in,'_',curr_session_block.file_name,'_min_num_',num2str(ii_file)];
            eval(['load ', current_filename]); % Load the current CSC data
            Current_file_timestamps = CSC.data.Timestamps_filtered_samples_saved;
            Current_file_timestamps_cleaned = Current_file_timestamps;
            
            %Load the associated artifacts file data
            current_artifact_filename = [filename_associated_artifact_file,'_',curr_session_block.file_name,'_min_num_',num2str(ii_file)];
            eval(['load ', current_artifact_filename]); % Load the current CSC data
           
           
    if (isempty(Current_file_timestamps)) % If there is an empty CSC file no need to clean it
        continue
    end
    
   
    for curr_chan_idx=1:4
        channel_cleaned{curr_chan_idx} = CSC.data.Samples_filtered_saved{used_channels(curr_chan_idx)};
    end
    
    
    for ii_channel = used_channels
                disp(['removing artifacts detected on channel #',num2str(ii_channel),' of file #', num2str(ii_file),' out of ',num2str(curr_session_block.file_num),' ', curr_session_block.session_name, ' files...'])
        
        for ii_artifact_on_this_channel = 1:size(Artifacts.data.Artifacts_start_end_IX,2)
            if (size(Artifacts.data.Artifacts_start_end_timestamps,1)>=ii_channel)
                current_artifact_start_end_IX = Artifacts.data.Artifacts_start_end_IX{ii_channel,ii_artifact_on_this_channel};
                if (length(current_artifact_start_end_IX)>1) % i.e., there is something there
                    current_artifact_start_timestamp = Current_file_timestamps(current_artifact_start_end_IX(1));
                    current_artifact_end_timestamp = Current_file_timestamps(current_artifact_start_end_IX(2));
                    idx_to_remove = find((Current_file_timestamps_cleaned>=current_artifact_start_timestamp)&(Current_file_timestamps_cleaned<=current_artifact_end_timestamp));
                    
                    % Now remove the samples & timestamps of those indexes from all the channels
                    for curr_chan_idx=1:4
                        channel_cleaned{curr_chan_idx}(idx_to_remove) = 0;
                    end
                    Current_file_timestamps_cleaned(idx_to_remove) = [];
                end
            end
        end
    end
    
    
   if p.remove_bad_frames==1
            % Remove the data for poor video data:
            if strcmp(curr_session_block.session_name,'behavioral')
                Current_file_timestamps_cleaned = Current_file_timestamps;
                for ii_frame = 1:length(bad_frames_start_end_timestamps)
                    current_frame_start_end = bad_frames_start_end_timestamps(ii_frame,:);
                    idx_to_remove = find((Current_file_timestamps_cleaned>=current_frame_start_end(1))&(Current_file_timestamps_cleaned<=current_frame_start_end(2)));
                    
                    % Now remove the samples & timestampls of those indexes from all the channels
                    for curr_chan_idx=1:4
                        channel_cleaned{curr_chan_idx}(idx_to_remove) = 0;
                    end
                    Current_file_timestamps_cleaned(idx_to_remove) = [];
                end
            end
            end
    
    % Detect crossing-threshold spikes and save their waveforms:
    Spike_waveforms_combined = [];
    Spike_timestamps_combined = [];
    
    % Initialize:
    for curr_chan_idx=1:4
        thres_cross_vec{curr_chan_idx} = zeros(1,length(channel_cleaned{curr_chan_idx}));
    end
    
    %             % look at what the spike threshold will cut off (in the first and last of minute of session) [Eyal]
    %             if (ii_file == 1 || ii_file == curr_session_block.file_num)						% process first and last file
    %                 threshold_fig = figure('name',['Threshold ',curr_session_block.session_name,' - Cleaned - F',num2str(ii_file)]);
    %                 for curr_chan_idx=1:4
    %                     subplot(4,1,curr_chan_idx)
    %                     plot(channel_cleaned{curr_chan_idx}(1:curr_session_block.plot_length)*HS_Conversion_Antenna(curr_chan_idx))
    %                     hold on
    %                     plot(1:curr_session_block.plot_length,ones(1,1:curr_session_block.plot_length)*Spike_threshold_uV_units,'r')
    %
    %                     if curr_chan_idx == 1
    %                         title(['Threshold ',curr_session_block.session_name,' - Cleaned - F', num2str(ii_file), ' - Converted to uV - Threshold = ', num2str(Spike_threshold_uV_units)])
    %                     end
    %                 end
    %
    %                 threshold_filename = [filename_cleaned_thresholds,curr_session_block.file_name,'_F', num2str(ii_file)];
    %                 saveas(threshold_fig, threshold_filename, 'fig');
    %                 saveas(threshold_fig, threshold_filename, 'jpg');
    %                 close(threshold_fig);
    %             end;
    %
    
    
    for curr_chan_idx= 1:1:size(p.active_channels,2)
        if curr_chan_idx == 2 %DBG
            flag = 1;
        end
        
        if p.active_channels(find(p.TT==ii_Tetrode),curr_chan_idx)==1 %i.e. - valid channel
            
            if p.include_negative_threshold == 1
                thres_cross_IX{curr_chan_idx} = find(channel_cleaned{curr_chan_idx}>Spike_threshold_uV_units | channel_cleaned{curr_chan_idx}<-Spike_threshold_uV_units);
            else
                thres_cross_IX{curr_chan_idx} = find(channel_cleaned{curr_chan_idx}>Spike_threshold_uV_units);
            end;
            
            % place '1' at threshold crossing:
            thres_cross_vec{curr_chan_idx}(thres_cross_IX{curr_chan_idx})=1; % place '1' at threshold crossing
            
            % Find the start and end segment of each spiking event on each channel of the tetrode:
            diff_thres_cross_vec{curr_chan_idx} = diff(thres_cross_vec{curr_chan_idx});
            SPK_start{curr_chan_idx} = find(diff_thres_cross_vec{curr_chan_idx} ==1)+1;
            SPK_End{curr_chan_idx} = find(diff_thres_cross_vec{curr_chan_idx} ==(-1))+1;
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
                counter = counter+ 1;
                temp_SPK_events = channel_cleaned{curr_chan_idx}(SPK_start{curr_chan_idx}(jj):SPK_End{curr_chan_idx}(jj));
                [temp_SPK_max,IX] = max(abs(temp_SPK_events));
                temp_IX_vec = SPK_start{curr_chan_idx}(jj):SPK_End{curr_chan_idx}(jj);
                max_IX = temp_IX_vec(IX);
                SPK_IX{curr_chan_idx}(1,jj) = max_IX; % Correct for the real timestamp relative to the entire recording session
            end
        else
            SPK_IX{curr_chan_idx} = [];
        end;
    end
    
    
    % Do coincedence detection to identify single spikes:
    % First find spikes which are too close apart to be considered as seperate spikes,
    % By two close I mean to have a temporal difference of < x_sep_spike_thres bins between them.
    SPK_All_Ch_combined = SPK_IX{1};
    
    for curr_chan_idx=2:4
        % First loop over the spikes detected on the first channel and compare their seperation from those detected on the second channel
        % Save only those which are seperated by a minimal number of bins to avoid counting the same spike twice
        for ii_spike = 1:length(SPK_All_Ch_combined)
            shared_IXs = [];
            current_spike_IX = SPK_All_Ch_combined(ii_spike);
            shared_IXs = find(abs(SPK_IX{curr_chan_idx} - current_spike_IX)<= x_sep_spike_thres);
            if ~isempty(shared_IXs) % i.e., the same spike is detected twice
                SPK_IX{curr_chan_idx}(shared_IXs) = [];
            else end
        end
        
        % Now that we do not have the same spikes on two channels we can merge the spike IXs of both channels,
        % as those represent different spikes
        SPK_All_Ch_combined = [SPK_All_Ch_combined,SPK_IX{curr_chan_idx}];
    end
    % sort them to preserve their temporal order:
    SPK_IXs_All_Ch_combined_sorted = sort(SPK_All_Ch_combined);
    
    
    %===============================================================
    % Extract the spike waveforms on all 4 channels of the tetrode:
    %===============================================================
    % Extract the waveform and timestamps of each spike (as in Neuralynx, the peak
    % will be the 8th sample out of over all 32 samples of
    %the spike shape vec.
    current_file_spike_counter = 0;
    for ii_spike = 1:length(SPK_IXs_All_Ch_combined_sorted)
        current_spike_max_IX = SPK_IXs_All_Ch_combined_sorted(ii_spike);
        if ((current_spike_max_IX+24<=length(Current_file_timestamps))& (current_spike_max_IX-7>0)) % i.e., we are NOT cutting the spike in the middle
            current_file_spike_counter = current_file_spike_counter + 1;
            SPK_timestamp(1,Last_Spike_IX+current_file_spike_counter) = Current_file_timestamps(current_spike_max_IX);
            %figure; %DBG
            %colors = 'brgk';
            %hold on
            for curr_chan_idx=1:4
                if p.active_channels(find(p.TT==ii_Tetrode),curr_chan_idx)==1 %i.e. - valid channel
                    SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter} = channel_cleaned{curr_chan_idx}(current_spike_max_IX-7:1:current_spike_max_IX+24);
                    %plot(SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter},colors(curr_chan_idx))
                    [temp_max(curr_chan_idx,1) temp_max(curr_chan_idx,2)] = max(SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter});
                else
                    SPK_waveforms{curr_chan_idx,Last_Spike_IX+current_file_spike_counter} =zeros(1,32);
                end;
            end
            
        end % if
    end % for - ii_spike
    Last_Spike_IX = Last_Spike_IX + current_file_spike_counter;
        end
clear CSC Artifacts ccc

end


%% Clean artifacts using the library of acceptable spike shapes - All spike waveforms:
%=====================================================================================

% Now we will throw away noisy theshold crossing using the library of
% acceptable spike shapes as reference
%  cd D:\Michael\Matlab ;
%
% % First we will need to normalize peak amplitude of each spike to a value
% % of '1' such that we can comapre it to the library of acceptable spike
% % shapes:
% Spike_waveforms_combined_and_norm = cell(size(SPK_waveforms,1),size(SPK_waveforms,2));
% for ii_spike = 1:size(SPK_waveforms,2)
%     for ii_channel = 1:size(SPK_waveforms,1)
%     temp_max = max(SPK_waveforms{ii_channel,ii_spike});
%     Spike_waveforms_combined_and_norm{ii_channel,ii_spike} = SPK_waveforms{ii_channel,ii_spike}./temp_max;
%     end
% end
count_gily = 1;
if p.do_comparison_acceptable_shapes==1 %if we  want to run comparison to librarty of acceptable shapes
    
    eval(['load ', library_file]); % Load the library of acceptable spike shapes;
    vector_of_accepted_spikes = zeros( 1, length(SPK_waveforms) ) + NaN ; % Initialize
    vector_of_max_r_values = zeros( 1, 1 ) + NaN ;
    
    for ii_spike = 1:size(SPK_waveforms,2)
        (ii_spike/size(SPK_waveforms,2))*100;
        spike_shape_4channels = zeros(32,4);
        for jj_channel = 1:size(SPK_waveforms,1)
            spike_shape_4channels(:,jj_channel) = SPK_waveforms{jj_channel,ii_spike}';
        end
        % Choose the channel # for which the spike has the largest height:
        %  [ stam  idx_channel_max_height ] = max( max( abs(spike_shape_4channels )) );
        [ stam  idx_channel_max_height ] = max( abs(spike_shape_4channels(8,:))) ;%GILY - changed to max of al channels just in the 8th spot (that was defined as peak)
        %[ stam  idx_channel_max_height ] = max( max( spike_shape_4channels ) ); %GILY - changed to abs in order to get reverse spike
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
        
      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Accept the spike shape ('1') if its correlation with ANY of the acceptable shapes
        % is >= r_threshold ; else, reject the spike ('0').
        
        % FOR DBG: plot waveform in 'r' if rejected and 'b' if accepted
        %     if (vector_of_sccepted_spikes( ii_spike )== 1)
        %     plot(spike_shape)
        %     else
        %     plot(spike_shape,'r')
        %     end
        %     title(['spike num = ',num2str(ii_spike), ' ,r val = ' num2str(vector_of_max_r_values( ii_spike ))])
        %     waitforbuttonpress
        
    end % End "Loop over spikes extracted from the Ntt file"
    % Find the IXs of the accepted spike waveforms and store the accepted and
    % not-accepted waveforms speratly.
    IX_accepted = find(vector_of_accepted_spikes ==1);
    IX_NO_accepted = find(vector_of_accepted_spikes ==0);
else % if we want to skip comparison to library
    vector_of_accepted_spikes=ones(1,length(SPK_waveforms));
    IX_accepted = 1:1:length(SPK_waveforms);
    IX_NO_accepted = [];
end;



% Extract the accepted (and not accepted) spikes and define new varialbes:
Spike_waveforms_accepted = cell(4,IX_accepted);
Spike_waveforms_NO_accepted = cell(4,IX_NO_accepted);

counter_accepted = 0;
counter_NO_accepted = 0;
for ii = 1:length(vector_of_accepted_spikes)
    if vector_of_accepted_spikes(ii) == 1;
        counter_accepted = counter_accepted + 1;
        
        for curr_chan_idx=1:4
            Spike_waveforms_accepted{curr_chan_idx,counter_accepted} = SPK_waveforms{curr_chan_idx,ii};
        end
    else
        counter_NO_accepted = counter_NO_accepted + 1;
        for curr_chan_idx=1:4
            Spike_waveforms_NO_accepted{curr_chan_idx,counter_NO_accepted} = SPK_waveforms{curr_chan_idx,ii};
        end
    end
end


%---------------------------------------------------------------
% Convert into NSE file for spike sorting - All spike waveforms
%---------------------------------------------------------------

%Note: We will NOT be using or saving the header at all

Timestamps_accepted_spikes = SPK_timestamp(IX_accepted);
%rotate the timestamps to be 1 x num_records
%timestamps = rot90(Timestamps_accepted_spikes);

% Clear un-needed variables
clear IX_NO_accepted IX_NO_accepted_sleep IX_accepted_sleep
clear SPK_All_Ch_combined
clear Spike_waveforms_NO_accepted
clear Spike_waveforms_NO_accepted_sleep Spike_waveforms_accepted_sleep Timestamps_accepted_spikes_sleep
clear VT VT_Parameters ccc spikes_sleep
clear vector_of_accepted_spikes vector_of_accepted_spikes_sleep

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

Timestamps_accepted_spikes_per_TT=Timestamps_accepted_spikes;
spikes_per_TT=spikes;

filename_out		= [dir_data_in, '\bat', num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode),'TEMP'];

save(filename_out, 'spikes_per_TT','Timestamps_accepted_spikes_per_TT')


end
