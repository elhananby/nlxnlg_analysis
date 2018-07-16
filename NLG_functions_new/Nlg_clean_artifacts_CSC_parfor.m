function Nlg_clean_artifacts_CSC_parfor (p)

%---------------------------
% Arseny Finkelstein (adapted from Michael Yartsev telemetry code) 03/04/2014
%--------------------------


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
pos_threshold=1500; % positive threshold for artifact removal (in uv)
neg_threshold=-1500; % positive threshold for artifact removal (in uv)

% Parameters extracted from the inclusion list (p) for each day

Day					= p.day;
Bat					= p.bat;

times_microsec_pre_sleep_session_1	= [p.S(1).start_time*1000 ;p.S(1).end_time*1000];
times_microsec_behav_session_2		= [p.S(2).start_time*1000 ;p.S(2).end_time*1000];
times_microsec_post_sleep_session_3 = [p.S(3).start_time*1000 ;p.S(3).end_time*1000];

% input
filename_in = [ p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\CSC_extracted_bat' num2str(Bat) '_Day' num2str(Day)]; % This is the prefix

% output
Output_Dir_Day=[p.path_dataout p.path_year_bat num2str(Day) '\'];
filename_out_dir = [Output_Dir_Day, 'Clean_artifacts\'];

% skip CSC Artifact cleaning if we had done this already
if isempty (dir(filename_out_dir))
    mkdir(filename_out_dir);
else
    disp('======================');
    disp (['Skipping CSC Artifact cleaning. Already extracted in ' filename_out_dir]);
    disp('======================');
    return;
end;


session_blocks(1).session_name = 'behavioral';
session_blocks(2).session_name = 'pre-sleep';
session_blocks(3).session_name = 'post-sleep';

session_blocks(1).file_name = 'behav';
session_blocks(2).file_name = 'pre_sleep';
session_blocks(3).file_name = 'post_sleep';

%load the number of CSC chunks to analyze
load ([p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\Num_of_analyzed_files.mat']);
session_blocks(1).file_num = Num_of_behav_analyzed_files;
session_blocks(2).file_num = Num_of_pre_sleep_analyzed_files;
session_blocks(3).file_num = Num_of_post_sleep_analyzed_files;

curr_session_block = session_blocks(1);

%-----------------------------------------
%% Find and extract artifacts from the session blocks for each tetrode
%-----------------------------------------
for ii_Tetrode=(p.TT)
    disp('======================');
    disp([' Cleaning CSC Artifact of Day #' num2str(Day) ' from Tetrode ' num2str(ii_Tetrode)]);
    disp('======================');
    
    filename_in = [ p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\CSC_extracted_bat' num2str(Bat) '_Day' num2str(Day) '_TT' num2str(ii_Tetrode),'_behav_min_num_']; % This is the prefix
    filename_out = [filename_out_dir, 'CSC_artifacts_bat', num2str(Bat), '_Day', num2str(Day) '_TT' num2str(ii_Tetrode),'_behav_min_num_']; % This is the prefix
    
     parfor ii_file = 1:curr_session_block.file_num % should be parfor
        Nlg_clean_artifacts_CSC_window(filename_in, filename_out, curr_session_block.file_num, ii_file)
        
     end
    %pre sleep
     filename_in = [ p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\CSC_extracted_bat' num2str(Bat) '_Day' num2str(Day) '_TT' num2str(ii_Tetrode),'_pre_sleep_min_num_']; % This is the prefix
    filename_out = [filename_out_dir, 'CSC_artifacts_bat', num2str(Bat), '_Day', num2str(Day) '_TT' num2str(ii_Tetrode),'_pre_sleep_min_num_']; % This is the prefix
    
     parfor ii_file = 1:curr_session_block.file_num % should be parfor
        Nlg_clean_artifacts_CSC_window(filename_in, filename_out, curr_session_block.file_num, ii_file)
        
    end
     %post sleep
     filename_in = [ p.path_datain p.path_year_bat 'CSC_extracted\'  num2str(Day) '\CSC_extracted_bat' num2str(Bat) '_Day' num2str(Day) '_TT' num2str(ii_Tetrode),'_post_sleep_min_num_']; % This is the prefix
    filename_out = [filename_out_dir, 'CSC_artifacts_bat', num2str(Bat), '_Day', num2str(Day) '_TT' num2str(ii_Tetrode),'_post_sleep_min_num_']; % This is the prefix
    
     parfor ii_file = 1:curr_session_block.file_num % should be parfor
        Nlg_clean_artifacts_CSC_window(filename_in, filename_out, curr_session_block.file_num, ii_file)
        
    end
    % print total  artifact samples
    % disp(['Total artifact time in ALL sessions recorded:', num2str(session_artifact_time)]);
    
end % end looping over Tetrodes

