function Nlg_Flight_detect_spikes_CSC_maya(p)

%% IMPORTANT to FIX: doesn't work in version newer than MATLAB 2012b because of the function cell

%---------------------------
% Maya G. - based on Arseny Finkelstein - 05/04/2014(adapted from Michael Yartsev telemetry code; Corrections - Eyal 23/6/2013)
%--------------------------


% Here we detect the spikes for the single tetrode on each antenna and save them.

%-----------------------------------------------------------------------
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
Bat					= p.bat;

used_channels=1:4;

% times_microsec_pre_sleep_session_1	= [p.S(1).start_time*1000 ;p.S(1).end_time*1000];
% times_microsec_behav_session_2		= [p.S(2).start_time*1000 ;p.S(2).end_time*1000];
% times_microsec_post_sleep_session_3 = [p.S(3).start_time*1000 ;p.S(3).end_time*1000];



%--------------------------------------------
% changed for every recording day.  for this batch file only

Spike_threshold_uV_units = p.Spike_threshold_uV_units;

%--------------------------------------------
% mostly constant.					for this batch file only

library_file = 'library_of_acceptable_spike_shapes.mat'; % File containing the Library Of Acceptable Spike Shapes


dir_data_in							= [p.path_dataout p.path_year_bat num2str(Day) '\'];
% %  filename_associated_VT_file			= [ p.path_dataout p.path_year_bat num2str(Day) '\' 'VT_flight_extracted_bat' num2str(Bat) '.mat']


% skip Spike extraction  if we had done this already
if ~isempty(dir(fullfile(dir_data_in, '*.NTT')))
    disp('======================');
    disp (['Skipping Spike detection from CSC. Already extracted in ' dir_data_in]);
    disp('======================');
    return;
end;

threshold_type = 1; % '1' - use a uV voltage threshold

%% Remove artifacts times and detect spikes
% Below we will ran, for each 1-min file, over the 4 channels of each
% antenna and remove all the CSC data during artifacts epochs from ALL
% channels of a single antenna.
% Then, we will run over the "artifact-Free" and extract the spikes on all
% 4-channels of a single antenna.

%-----------------------------------------

%% Extract data from the session blocks

% % session_blocks(1).session_name = 'pre-sleep';
% % session_blocks(2).session_name = 'behavioral';
% % session_blocks(3).session_name = 'post-sleep';
% %
% % session_blocks(1).file_name = 'pre_sleep';
% % session_blocks(2).file_name = 'behav';
% % session_blocks(3).file_name = 'post_sleep';

%load the number of CSC chunks to analyze
load ([dir_data_in '\CSC_extracted\' , 'Num_of_analyzed_files.mat']);
session_blocks(1).file_num = Num_of_analyzed_files;
curr_session_block = session_blocks(1);



% % session_blocks(1).file_num = Num_of_pre_sleep_analyzed_files;
% % session_blocks(2).file_num = Num_of_behav_analyzed_files;
% % session_blocks(3).file_num = Num_of_post_sleep_analyzed_files;
% %
% %
% % session_blocks(1).plot_length = 10000;
% % session_blocks(2).plot_length = 50000;
% % session_blocks(3).plot_length = 10000;

% % eval(['load ', filename_associated_VT_file]);
% % bad_frames_start_end_timestamps = VT.bad_frames_start_end_timestamps;


Timestamps_accepted_spikes_TT=[];% Initialize
spikes_TT=[];% Initialize

tetrodes=p.use_tetrodes;
TT_counter=0;

% Loop over tetrodes, a TEMP file is saved for each tetrode and loaded in
% the next step:
parfor ii_Tetrode=(p.use_tetrodes)
    
    disp('======================');
    disp([' Detecting Spikes from CSC. #' num2str(Day) ' from Tetrode ' num2str(ii_Tetrode)]);
    disp('======================');
    
    if ~isempty(dir(fullfile(dir_data_in ,sprintf('*bat*TT%d*TEMP.mat',ii_Tetrode))))
        disp('======================');
        disp (['Skipping 1st phase of Spike detection from CSC. Already exists in ' dir_data_in]);
        disp('======================');
    else
        
        Nlg_Flight_detect_spikes_CSC_per_TT(p, ii_Tetrode)
    end
    
end; %TTs

for ii_Tetrode=(p.use_tetrodes)
    filename_out		= [dir_data_in, '\bat', num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode),'TEMP.mat'];
    load(filename_out);
    Timestamps_accepted_spikes_TT{ii_Tetrode} = Timestamps_accepted_spikes_per_TT;
    spikes_TT{ii_Tetrode}=spikes_per_TT;
    delete(filename_out);
    clear Timestamps_accepted_spikes_per_TT spikes_per_TT
end


%%  Read Timestamp data and use Coincidence-Detection across Tetrodes to eliminate artifacts (just as we do in the wired case): --------

idx_coincidence_vec{tetrodes(end)} = [];  % Initialize this variable (for later)

% Find coincidence-detection events = coincidence-detection on a millisecond-scale:
for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{tetrodes(1)}), % Loop over the spikes of the FIRST tetrode
    t_spike = Timestamps_accepted_spikes_TT{tetrodes(1)}(ii_spikes);
    idx{tetrodes(1)} = ii_spikes; % Index of this spike
    
    for ii_file = 2:1:length(tetrodes) % Loop over the other tetrodes, finding the coincidnce-detection
        idx{tetrodes(ii_file)} = find( abs( Timestamps_accepted_spikes_TT{tetrodes(ii_file)} - t_spike ) <= coincidence_window ); % THE COINCIDENCE DETECTION
    end
    
    % Check that the coincidence-detection occurred on ALL tetrodes:
    test_variable = 1;
    for ii_file =  tetrodes % Loop over all tetrode files
        if ( isempty( idx{ii_file} ) ), % If there is NO coincidence-detection
            test_variable = 0 ;
        end
    end
    
    if length(tetrodes)>=3 %If there are at least 3 tetrodes
        % Check that the temporal separation between the OTHER two tetrodes also meets the coincidence-detection criterion
        % (this is needed since the spikes on the OTHER two tetrodes may be up to twice-the-time apart, in principle!!! ):
        if ( test_variable == 1 ), % If we passed the first test
            if ( abs( Timestamps_accepted_spikes_TT{tetrodes(2)}(idx{tetrodes(2)}(1)) - Timestamps_accepted_spikes_TT{tetrodes(3)}(idx{tetrodes(3)}(1)) ) > coincidence_window ),
                test_variable = 0 ; % Reset test_variable if the time between spikes on the OTHER two tetrodes is too large
            end
        end
    end;
    
    % Save the indexes of the coincidence-detection "spikes" to be REMOVED:
    if ( test_variable == 1 ), % If all indexes exist = there IS a coincidence detection on ALL tetrodes
        for ii_file = tetrodes % Loop over all tetrode files
            idx_coincidence_vec{ii_file} = [ idx_coincidence_vec{ii_file}  idx{ii_file} ];
        end
    end
    if mod(ii_spikes,1000) == 0
        disp(['Processing coindicedence detection across TT- ' num2str((ii_spikes/length(Timestamps_accepted_spikes_TT{tetrodes(1)}))*100)])
    end
    
end %end looping over spikes for coincidence detection


%% Save the data in Neuralynx NTT files for three cases:
% (1) All the data.
% (2) Sleep sessions only (pre and post-sleep sessions).
% (3) Behav session only.

for ii_Tetrode= (p.use_tetrodes)
    
    filename_out		= [dir_data_in, '\bat', num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode)];
    filename_out_NLG	= [dir_data_in, '\bat', num2str(Bat), '_Day', num2str(Day), '_TT', num2str(ii_Tetrode),'_NLG'];

    % load generic header to add to NTT files
    load('ntt_generic_header.mat');
    
    %removing spikes because of coincidence detection
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{ii_Tetrode};
    Timestamps_accepted_spikes(idx_coincidence_vec{ii_Tetrode})=[];
    spikes=spikes_TT{ii_Tetrode};
    spikes(:,:,idx_coincidence_vec{ii_Tetrode})=[];
    %only export timestamps and data points - Full session
    FieldSelection = [1 0 0 0 1 1];
    
    % Maya - transform spikes to NLX clk - 
    Timestamps_accepted_spikes_NLX = polyval(p.polyfit_Nlg2Nlx_microsec{p.nlg_part}, Timestamps_accepted_spikes);
    Mat2NlxSpike( [filename_out,'.NTT'], 0, 1, 1, FieldSelection, Timestamps_accepted_spikes_NLX , spikes, Header);
    Mat2NlxSpike( [filename_out_NLG,'.NTT'], 0, 1, 1, FieldSelection, Timestamps_accepted_spikes , spikes, Header);

    
    %     %only export timestamps and data points - Sleep ONLY
    %     Sleep_spike_IXs_pre_sleep = find((Timestamps_accepted_spikes>times_microsec_pre_sleep_session_1(1))&(Timestamps_accepted_spikes<times_microsec_pre_sleep_session_1(2)));
    %     if~isempty(times_microsec_post_sleep_session_3)
    %         Sleep_spike_IXs_post_sleep = find((Timestamps_accepted_spikes>times_microsec_post_sleep_session_3(1))&(Timestamps_accepted_spikes<times_microsec_post_sleep_session_3(2)));
    %     else
    %         Sleep_spike_IXs_post_sleep = [];
    %     end
    %     Sleep_spike_IXs = [Sleep_spike_IXs_pre_sleep,Sleep_spike_IXs_post_sleep];
    %     Timestamps_accepted_spikes_sleep = Timestamps_accepted_spikes(Sleep_spike_IXs);
    %     spikes_sleep = spikes(:,:,Sleep_spike_IXs);
    %     FieldSelection = [1 0 0 0 1 1];
    %     %Mat2NlxTT( [filename_out_sleep,'.NTT'], 0, 1, 1, length(Timestamps_accepted_spikes_sleep), FieldSelection, Timestamps_accepted_spikes_sleep , spikes_sleep, Header);
    %     Mat2NlxSpike( [filename_out_sleep,'.NTT'], 0, 1, 1, FieldSelection, Timestamps_accepted_spikes_sleep , spikes_sleep, Header);
    %
    %     %only export timestamps and data points - Behav ONLY
    %     Behav_spike_IXs = find((Timestamps_accepted_spikes>times_microsec_behav_session_2(1))&(Timestamps_accepted_spikes<times_microsec_behav_session_2(2)));
    %     Timestamps_accepted_spikes_behav = Timestamps_accepted_spikes(Behav_spike_IXs);
    %     spikes_behav = spikes(:,:,Behav_spike_IXs);
    %     FieldSelection = [1 0 0 0 1 1];
    %     %Mat2NlxTT( [filename_out_behav,'.NTT'], 0, 1, 1, length(Timestamps_accepted_spikes_behav), FieldSelection, Timestamps_accepted_spikes_behav , spikes_behav, Header);
    %     Mat2NlxSpike( [filename_out_behav,'.NTT'], 0, 1, 1, FieldSelection, Timestamps_accepted_spikes_behav , spikes_behav, Header);
    %
    
    
end % end looping over tetrodes

