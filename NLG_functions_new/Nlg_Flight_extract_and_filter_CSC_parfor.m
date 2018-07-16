function  Nlg_Flight_extract_and_filter_CSC_parfor (p)

%-----------------------------------------------------------------------
%% PARAMETERS:
Num_of_pre_sleep_analyzed_files=[]; %initializing
Num_of_behav_analyzed_files=[]; %initializing
Num_of_post_sleep_analyzed_files=[]; %initializing
Window_size_used = 2*60*10^6; % we will run over the data in 2-min windows but save
%only the central 1-min to avoid edge effects in the filtering process (see below)
Window_size_saved = 60*10^6; % we will run over the data in 1-min windows
Filter__N_order = 300 ; %  Same as we use for filtering ripples
Filter_FreqsHz = [600 6000]; % Filter for spikes
amplitude_divide_factor = 1;
%amplitude_divide_factor = 20; % This is to bring the values closer to the sin fitting model in order to aid faster
% assigment of the variables.

% Parameters extracted from the inclusion list (p) for each day

Day					= p.day;
Bat					= p.bat;

Output_Dir_Root		= p.path_dataout;

times_microsec_pre_sleep_session_1	= [p.S(1).start_time ;p.S(1).end_time];
times_microsec_behav_session_2		= [p.S(2).start_time ;p.S(2).end_time];
times_microsec_post_sleep_session_3 = [p.S(3).start_time ;p.S(3).end_time];

% input
Input_Dir_Recording = [p.path_datain p.path_year_bat p.path_day_dir '\'];

Recording_directory_general = [Input_Dir_Recording 'nlx\'];
Recording_directory_pre_sleep_session	= Recording_directory_general;	% We can seperate the recording directories for the cases where we had to re-start the telemetry computer during recordings
Recording_directory_behav_session		= Recording_directory_general;
Recording_directory_post_sleep_session	= Recording_directory_general;

t_throw_away_data	= [];	% Throw away times when the bat DID NOT BEHAVE (see: hippo_clean_artifacts_batch_Michael.m)

% output
csc_dir					 = [ p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\'];
file_dir_out_for_figure  = [csc_dir];

% skip CSC extraction and exit function if we have already extracted the data
if isempty (dir(csc_dir)) 
    mkdir(csc_dir);
else
    disp('======================');
    disp (['Skipping CSC extraction. Already extracted in ' csc_dir]);
    disp('======================');
    return;
end;




non_active_channels = [];

%% Extract data for each Tetrode

for ii_Tetrode=(p.TT)
    
    disp('======================');
    disp([' Extracting CSC of Day #' num2str(Day) ' from Tetrode ' num2str(ii_Tetrode)]);
    disp('======================');
    
    filename_out_save = [csc_dir, 'CSC_extracted_bat', num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode)];
    
    Active_channels = 1:4;
    
    %-----------------------------------------
    %% Extract data from the behavioral session
    %-----------------------------------------
    
    
    num_windows = floor(diff(times_microsec_behav_session_2)/Window_size_saved)-1; % Note that we loss ~ 0-1 min
    % of data on this, but I do this because we are taking a window for the
    % filtering twice the size of the window we wish to save and hence if we
    % were to use the maximal number of windown (1 more then what we currently
    % use) we would be using data outside this session.
    %num_windows = floor(diff(times_microsec_behav_session_2)/Window_size_saved);
    
    Num_of_behav_analyzed_files=num_windows;
    
   parfor ii_window = 1:num_windows
        filename_out{ii_window}= [filename_out_save,'_behav_min_num_',num2str(ii_window),'.mat'];
        NLG_filter_CSC_window(filename_out{ii_window}, times_microsec_behav_session_2, ii_Tetrode, ii_window, Recording_directory_general, num_windows)
    end
    
    
    %% Extract data from the pre-sleep session
    %-----------------------------------------
    
    if  ~isempty(times_microsec_pre_sleep_session_1) % if this session was even was ran

        num_windows = floor(diff(times_microsec_pre_sleep_session_1)/Window_size_saved)-1;
        Num_of_pre_sleep_analyzed_files=num_windows;

        parfor ii_window = 1:num_windows
                   filename_out{ii_window} = [filename_out_save,'_pre_sleep_min_num_',num2str(ii_window),'.mat'];

                    NLG_filter_CSC_window(filename_out{ii_window}, times_microsec_pre_sleep_session_1, ii_Tetrode, ii_window, Recording_directory_general, num_windows)

                  end
        
    else end
    
    
    %% Extract data from the post-sleep session
    %-----------------------------------------
    
    if  ~isempty(times_microsec_post_sleep_session_3) % if this session was even was ran
        num_windows = floor(diff(times_microsec_post_sleep_session_3)/Window_size_saved)-1;
        Num_of_post_sleep_analyzed_files=num_windows;
        parfor ii_window = 1:num_windows
                        filename_out{ii_window} = [filename_out_save,'_post_sleep_min_num_',num2str(ii_window),'.mat'];

                    NLG_filter_CSC_window(filename_out{ii_window}, times_microsec_post_sleep_session_3, ii_Tetrode, ii_window, Recording_directory_general, num_windows)
        end
        
    else end
    
end % end looping over Tetrodes


Num_of_analyzed_files.Num_of_pre_sleep_analyzed_files = Num_of_pre_sleep_analyzed_files;
Num_of_analyzed_files.Num_of_behav_analyzed_files = Num_of_behav_analyzed_files;
Num_of_analyzed_files.Num_of_post_sleep_analyzed_files = Num_of_post_sleep_analyzed_files;

save ([csc_dir 'Num_of_analyzed_files'], '-struct' ,'Num_of_analyzed_files');