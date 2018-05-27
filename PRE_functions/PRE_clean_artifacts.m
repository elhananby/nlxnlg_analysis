function PRE_clean_artifacts(p,datadir_out)

% Cleaning artifacts in Neuralynx Ntt files; this script is run via:  hippo_clean_artifacts_batch_Michael.m .
% The cleaning is based on:
% (1) Coincidence Detection on a sub-millisecond scale between ALL tetrodes (eliminates artefacts).
% (2) Session times (eliminate the inter-session times).
% (3) Comparison of recorded waveforms to a Library of Acceptable Spikes.


%-----------------
% Dori Derdikman, adapted from Michael Y. and Maya G.
%-----------------

% start batch processing

datadir_in = fullfile(p.path_datain, p.data_dir);
TT_vec = p.use_tetrodes;

filenames_in = {};
for TT = TT_vec
    filenames_in{end+1}=sprintf('TT%d.Ntt',TT);
end
files_to_export = (ismember(TT_vec,p.TT));
%Meaning, this vector tells us which TT (from the variable "filenames_in" we want to analyze.

% if ntt files were already created after clean_artifacts, do not perform
% function

files_interesting = find( files_to_export == 1 ) ; % These are the files to be processed and exported
prefix_for_saving = sprintf('animal%d_Day%d_%d_', p.animal, p.day, p.experiment);

% check if all files already exist
if exist(datadir_out, 'dir')
    nttFiles = subdir(sprintf('%s\\%s*.ntt', datadir_out, prefix_for_saving));
    if ~isempty(nttFiles) && length(nttFiles) == sum(files_to_export)
        fprintf('Skipping spike file creation, %s already exists\n', datadir_out);
        return;
    end
end

t_throw_away_data = p.throw_away_times*1e6; % Throw away times when the mouse DID NOT BEHAVE.

% ======= Parameters: ==================
% Parameters are defined in the batch-script:  hippo_clean_artifacts_batch.m
% ======================================

% ------------  Some parameters that do NOT change between different data-files: -------------

coincidence_window = 500 ; % Duration of coincidence-detection window (in microsec)

library_file = ... % File containing the Library Of Acceptable Spike Shapes
    fullfile(p.path_dataout,'library_of_acceptable_spike_shapes.mat');

% Reject recorded waveforms that have correlation < r_threshold with ALL of the "acceptable spike shapes"
r_threshold = p.r_threshold;

disp(' ');
disp(['Processing ', datadir_in]);
disp(' ');

% --------- STEP 1:  Take only spikes that occured within sessions: ---------------

files_interesting = find( files_to_export == 1 ) ; % These are the files to be processed and exported

for ii_file = files_interesting % Loop over the "interesting" files only
    
    filename_full_IN = fullfile(datadir_in, filenames_in{ii_file});
    %[datadir_in, '\', filenames_in{ii_file}];
    
    FieldSelection = [1 0 0 0 0] ; % Only get Timestamps
    ExtractionMode = 1 ; % Will read all the records
    ExtractionModeArray = [] ; % Will read all the records
    [Timestamps, NlxHeader] = ...
        Nlx2MatSpike( filename_full_IN, FieldSelection, 1, ExtractionMode, ExtractionModeArray );
    
    TotalSpikeNumber{ii_file} = length(Timestamps);
    idx_sessions = []; % Initialize
    
    start_times = [p.S.start_time];
    end_times = [p.S.end_time];
    nsessions = length(start_times);
    
    for i=1:nsessions
        idx_sessions = [idx_sessions ...
            find( Timestamps > start_times(i) & Timestamps < end_times(i) )];
    end
    
    idx_sessions = unique( idx_sessions ); % This is a Sorted vector
    
    % If I need to throw away data WITHIN a session:
    % Why is this done while looping over tetrodes? It's a global thing
    [number_of_throw_away_events, throw_away_times_of_each_event] = size(t_throw_away_data);
    
    idx_remove_from_the_idx_sessions_variable = [];% Initialize
    
    if ( throw_away_times_of_each_event == 2 )
        for throw_away_event_num = 1:number_of_throw_away_events
            
            throw_start_idx = t_throw_away_data(throw_away_event_num, 1);
            throw_end_idx = t_throw_away_data(throw_away_event_num, 2);
            
            idx_remove_from_the_idx_sessions_variable = [idx_remove_from_the_idx_sessions_variable,...
                find( Timestamps( idx_sessions ) > throw_start_idx & ...
                Timestamps( idx_sessions ) < throw_end_idx)];
        end
        
        idx_sessions( idx_remove_from_the_idx_sessions_variable ) = [];
    end
    
    idx_inter_sessions_all{ii_file} = 1:TotalSpikeNumber{ii_file};
    idx_inter_sessions_all{ii_file}(idx_sessions) = []; %
    % This variable will now store the spikes to REMOVE = the spikes that were BETWEEN sessions
    
end


% ------- STEP 2:  Read Timestamp data and use Coincidence-Detection to eliminate artifacts: --------

for ii_file = 1:length(filenames_in) % Loop over the tetrode SC files (Ntt files): Read Timestamps
    
    filename_full_IN = fullfile(datadir_in, filenames_in{ii_file});
    
    FieldSelection = [1 0 0 0 0] ; % Only get Timestamps
    ExtractionMode = 1 ; % Will read all the records
    ExtractionModeArray = [] ; % Will read all the records
    [Timestamps, NlxHeader] = ...
        Nlx2MatSpike( filename_full_IN, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
    
    Timestamps_all{ii_file} = Timestamps; % timestamps in microsec (holds timing of spikes for each tetrode)
    
    idx_coincidence_vec{ii_file} = [];  % Initialize this variable (for later)
    
end

% Find coincidence-detection events = coincidence-detection on a millisecond-scale:

for ii_spikes = 1:length(Timestamps_all{1}) % Loop over the spikes of the FIRST tetrode
    t_spike = Timestamps_all{1}(ii_spikes);
    idx{1} = ii_spikes; % Index of this spike
    for ii_file = 2:length(filenames_in) % Loop over the other tetrodes, finding the coincidnce-detection
        idx{ii_file} = find( abs( Timestamps_all{ii_file} - t_spike ) <= coincidence_window ); % THE COINCIDENCE DETECTION
    end
    
    % Check that the coincidence-detection occurred on ALL tetrodes:
    test_variable = 1;
    for ii_file = 1:length(filenames_in) % Loop over all tetrode files
        if ( isempty( idx{ii_file} ) ) % If there is NO coincidence-detection
            test_variable = 0 ;
        end
    end
    
    % Check that the temporal separation between the OTHER two tetrodes also meets the coincidence-detection criterion
    % (this is needed since the spikes on the OTHER two tetrodes may be up to twice-the-time apart, in principle!!! ):
    if ( test_variable == 1 ) % If we passed the first test
        if length(Timestamps_all) > 2 && ...
                ( abs( Timestamps_all{2}(idx{2}(1)) - Timestamps_all{3}(idx{3}(1)) ) > coincidence_window )
            test_variable = 0 ; % Reset test_variable if the time between spikes on the OTHER two tetrodes is too large
        end
    end
    
    % Save the indexes of the coincidence-detection "spikes" to be REMOVED:
    if ( test_variable == 1 ) % If all indexes exist = there IS a coincidence detection on ALL tetrodes
        for ii_file = 1:length(filenames_in) % Loop over all tetrode files
            idx_coincidence_vec{ii_file} = [ idx_coincidence_vec{ii_file}  idx{ii_file} ];
        end
    end
    
end

clear  Samples ChanNum CellNumbersSpikeSorting NumValidSamples idx_sessions Timestamps_all idx_spikes_remove vector_of_sccepted_spikes vector_of_max_r_values ; % Free up Memory



% --------- STEP 3:  Remove waveforms that do NOT match the Library Of Acceptable Spikes: ---------

for ii_file = files_interesting % Loop over the "interesting" files only
    
    idx_not_acceptable_shape{ii_file} = [];
    filename_full_IN = fullfile(datadir_in, filenames_in{ii_file});
    
    FieldSelection = [0 0 0 0 1]; % Will read ONLY the Samples (waveforms)
    ExtractionMode = 2; % Extract Record Index Range. (Count starts at '0')
    ExtractionModeArray = [1 TotalSpikeNumber{ii_file}];
    
    Samples = Nlx2MatSpike( filename_full_IN, FieldSelection, 0, ExtractionMode, ExtractionModeArray ) ;
    
    % Comparison to the Library Of Acceptable Spikes:
    [ vector_of_accepted_spikes,  vector_of_max_r_values_tmp ] = ...
        library_of_acceptable_spike_shapes__run_algorithm(ii_file, library_file, Samples, r_threshold);
    
    idx_not_acceptable_shape{ii_file} = find( vector_of_accepted_spikes == 0 );
    clear Samples
end

% ------ STEPS 1+2+3:  Combine spikes-to-be-removed from Steps 1+2+3, remove them, and SAVE: ------

TT_index = 0;

for ii_file = files_interesting % Loop over the "interesting" files only
    
    TT_index = TT_index+1;
    TT = p.TT(TT_index);
    idx_spikes_remove{ii_file} = ...  % Spikes to remove
        unique( [...
        idx_coincidence_vec{ii_file} , ... % Coincidence-detection spikes (coincidence-detection on all Tetrodes)
        idx_inter_sessions_all{ii_file},  ... % Inter-session times (throw those away)
        idx_not_acceptable_shape{ii_file} ... % Spikes with waveforms that do NOT match the Library Of Acceptable Spikes
        ] );
    
    idx_spikes_choose{ii_file} = 1 : TotalSpikeNumber{ii_file};  % The complementary list = Spikes to Choose
    idx_spikes_choose{ii_file}(idx_spikes_remove{ii_file}) = [];
    
    % Now do the actual cleaning of the data -- by reading the Ntt file again, including the Samples
    % variable, but taking only the desired ("non-cleaned") records:
    
    %% read original file
    filename_full_IN = fullfile(datadir_in, filenames_in{ii_file});
    
    filename_OUT = sprintf('%s_TT%i.ntt' ,prefix_for_saving, TT);
    
    filename_full_OUT = fullfile(datadir_out, filename_OUT);
    
    if ~exist(fileparts(filename_full_OUT), 'dir')
        mkdir(fileparts(filename_full_OUT));
    end
    
    FieldSelection = [1 1 1 1 1];
    ExtractionMode = 3;
    [Timestamps, ScNumbers, CellNumbersSpikeSorting, Features, Samples, NlxHeader] =...
        Nlx2MatSpike( filename_full_IN, FieldSelection, 1, ExtractionMode, idx_spikes_choose{ii_file});
    
    % Save the cleaned data into a new Ntt file (save ALL the variables):   
    features_for_header = load('Header_Features.mat');
    fet_for_head = strtrim(features_for_header.Header2);
    
    NlxHeader(end-7:end) = [];
    
    NlxHeader = [NlxHeader; fet_for_head];
    
    Mat2NlxSpike( filename_full_OUT, 0, 1, 1,...
        [1 1 1 1 1 1], Timestamps, ScNumbers, CellNumbersSpikeSorting,...
        Features, Samples, NlxHeader);
    
    filename_PARAMS_out = fullfile(datadir_out, '\params_clean_artifacts.mat');
    
    load('D:\Scripts\nlx_analysis\library_of_acceptable_spike_shapes.mat')
    
    functional_tetrode_channels = p.active_channels(TT_index,:); % '0' = dys-functional channel, '1' otherwise (this variable will be stored, but will NOT affect the output Ntt file)
    save(filename_PARAMS_out, 'coincidence_window' , 'r_threshold' , 'functional_tetrode_channels' , 'library_of_acceptable_spike_shapes')
    
    
end
end
