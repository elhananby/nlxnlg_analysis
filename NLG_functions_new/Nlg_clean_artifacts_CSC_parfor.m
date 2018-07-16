function Nlg_clean_artifacts_CSC_parfor (p)
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
pos_threshold = 1500; % positive threshold for artifact removal (in uv)
neg_threshold = -1500; % positive threshold for artifact removal (in uv)

% Parameters extracted from the inclusion list (p) for each day

times_microsec_behav_session_2		= [p.nlg.start_time*1000 p.nlg.end_time*1000];

Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');	% We can seperate the recording directories for the cases where we had to re-start the telemetry computer during recordings

csc_dir	= fullfile(Recording_directory_general, 'CSC_extracted\');


% output
Output_Dir_Day = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');
filename_out_dir = [Output_Dir_Day, '\Clean_artifacts\'];

% skip CSC Artifact cleaning if we had done this already
if isempty (dir(filename_out_dir))
    mkdir(filename_out_dir);
% else
%     disp('======================');
%     disp (['Skipping CSC Artifact cleaning. Already extracted in ' filename_out_dir]);
%     disp('======================');
%     return;
end


session_blocks(1).session_name = 'behavioral';

session_blocks(1).file_name = 'behav';

%load the number of CSC chunks to analyze
load ([csc_dir '\Num_of_analyzed_files.mat']);
session_blocks(1).file_num = Num_of_behav_analyzed_files;

curr_session_block = session_blocks(1);

%-----------------------------------------
%% Find and extract artifacts from the session blocks for each tetrode
%-----------------------------------------
for ii_Tetrode=(p.TT)
    fprintf('Processing Tetrode %i\n', ii_Tetrode);
    filename_in = [csc_dir, 'CSC_extracted_animal', num2str(p.animal), '_Day', num2str(p.day), '_TT', num2str(ii_Tetrode), '_behav_min_num_'];
    filename_out = [filename_out_dir, 'CSC_artifacts_animal', num2str(p.animal), '_Day', num2str(p.day) '_TT' num2str(ii_Tetrode),'_behav_min_num_']; % This is the prefix
    
     for ii_file = 1:curr_session_block.file_num % should be parfor
        Nlg_clean_artifacts_CSC_window(filename_in, filename_out, curr_session_block.file_num, ii_file)
     end

    
end % end looping over Tetrodes

