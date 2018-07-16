function NLG_filter_CSC_window(filename_out, times_microsec_behav_session_2, ii_Tetrode, ii_window, Recording_directory_general, num_windows)

%% PARAMETERS:
Num_of_pre_sleep_analyzed_files=[]; %initializing
Num_of_behav_analyzed_files=[]; %initializing
Num_of_post_sleep_analyzed_files=[]; %initializing
Window_size_used = 2*60*10^6; % we will run over the data in 2-min windows but save
%only the central 1-min to avoid edge effects in the filtering process (see below)
Window_size_saved = 60*10^6; % we will run over the data in 1-min windows
Filter__N_order = 300 ; %  Same as we use for filtering ripples
Filter_FreqsHz = [600 6000]; % Filter for spikes
CSC_Sampling_Rate_Hz = 28409.09;%of the NEUROLOGGER
amplitude_divide_factor = 1;
%amplitude_divide_factor = 20; % This is to bring the values closer to the sin fitting model in order to aid faster
% assigment of the variables.
Active_channels = 1:4;
non_active_channels = [];

t_throw_away_data	= [];	% Throw away times when the bat DID NOT BEHAVE (see: hippo_clean_artifacts_batch_Michael.m)

disp(['Behavioral data of minute No. ',num2str(ii_window),' out of ',num2str(num_windows),' are now being processed...'])
Tstart = times_microsec_behav_session_2(1) + Window_size_saved*(ii_window-1);
Tend = Tstart + Window_size_used;

for ii_channel = Active_channels
    
    channel_to_extract= (ii_channel-1) + 4*(ii_Tetrode-1);
    
    filename_CSC_EEG_in = fullfile(Recording_directory_general,['CSC',num2str(channel_to_extract),'.ncs']);
    %disp(['Behavioral data of channel ',num2str(ii_channel),' are now being processed...'])
    
    Samples_filtered = [];
    % Cut out the current samples
    [Samples_current, Timestamps_current] = Nlg_CutSamples_CSC(filename_CSC_EEG_in, Tstart, Tend);
    
    %Filter the current samples
    Wn = Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
    win_filter = fir1( Filter__N_order, Wn, 'bandpass' ); % Filter parameters are defined above
    
    if length(Samples_current) < 3*length(win_filter)
        warning('Not enough data in minute %d',ii_window)
        
        CSC.params.Window_size_used = Window_size_used ;
        CSC.params.Window_size_saved = Window_size_saved ;
        CSC.params.Filter__N_order = Filter__N_order ;
        CSC.params.Filter_FreqsHz = Filter_FreqsHz ;
        CSC.params.CSC_Sampling_Rate_Hz = CSC_Sampling_Rate_Hz ;
        CSC.params.Recording_directory = Recording_directory_general;
        CSC.params.non_active_channels = non_active_channels; % These will get a value of zero
        CSC.data.Samples_filtered_saved = [] ;
        CSC.data.Timestamps_filtered_samples_saved = [] ;
        
        %CSC.data.Samples_original_unfiltered = Samples_current ;
        %         CSC.data.HS_conversion = HS_conversion;
        
        save(filename_out, 'CSC')

        return
    end
    
    Samples_filt = filtfilt( win_filter, 1, Samples_current ); % Use UN-DECIMATED data

    % Find the IXs that we wrote to ourselfs (during the recording itself) to
    % throw away from the behav session:
%     if ~isempty(t_throw_away_data)
%         IXs_to_throw_away = find((Timestamps_current>t_throw_away_data(1))&((Timestamps_current<t_throw_away_data(2))));
%         if ~isempty(IXs_to_throw_away)
%             Samples_filt(IXs_to_throw_away) = [];
%             Timestamps_current(IXs_to_throw_away) = [];
%         end
%     end
    
    
    if length(Timestamps_current)>3 % i.e., we didn't throw away all the data for this window in the previous stage
        % Save the data only for the central portion of the filtered
        % segment (to avoid edge effects of the filtering):
        Samples_filtered_saved{ii_channel} = Samples_filt(length(Timestamps_current)*0.25:length(Timestamps_current)*0.75 - 1); % Save only the central 1-min of the samples (to avoid edge effects of the filtering)
        
        % The timestamps below are of ocurse the same for all channels so I
        % will save it only once:
        Timestamps_saved = Timestamps_current(length(Timestamps_current)*0.25:length(Timestamps_current)*0.75 - 1); % Save only the central 1-min of the samples (to avoid edge effects of the filtering)
        
    else
        Samples_filtered_saved{ii_channel} = [];
        Timestamps_saved = [];
    end
    
    % If this channel is in-active, simply place zero's in the samples
    % variable.
    if ismember(ii_channel,non_active_channels)
        Samples_filtered_saved{ii_channel}= zeros(1,length(Samples_filtered_saved{ii_channel}));
    else end
     
end

CSC.params.Window_size_used = Window_size_used ;
CSC.params.Window_size_saved = Window_size_saved ;
CSC.params.Filter__N_order = Filter__N_order ;
CSC.params.Filter_FreqsHz = Filter_FreqsHz ;
CSC.params.CSC_Sampling_Rate_Hz = CSC_Sampling_Rate_Hz ;
CSC.params.Recording_directory = Recording_directory_general;
CSC.params.non_active_channels = non_active_channels; % These will get a value of zero
CSC.data.Samples_filtered_saved = Samples_filtered_saved ;
CSC.data.Timestamps_filtered_samples_saved = Timestamps_saved ;

%CSC.data.Samples_original_unfiltered = Samples_current ;
%         CSC.data.HS_conversion = HS_conversion;

save(filename_out, 'CSC')
clear CSC Samples_filt Samples_current Timestamps_current