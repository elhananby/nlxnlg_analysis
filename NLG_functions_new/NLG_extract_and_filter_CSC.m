function  NLG_extract_and_filter_CSC(p)
%% PARAMETERS:
%only the central 1-min to avoid edge effects in the filtering process (see below)
Window_size_saved = 60*10^6; % we will run over the data in 1-min windows

times_microsec_behav_session_2 = [p.nlg.start_time ;p.nlg.end_time];

% input
Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');	% We can seperate the recording directories for the cases where we had to re-start the telemetry computer during recordings

% output
csc_dir	= fullfile(Recording_directory_general, 'CSC_extracted\');

% skip CSC extraction and exit function if we have already extracted the data
if isempty (dir(csc_dir))
    mkdir(csc_dir);
else
    disp('======================');
    disp (['Skipping CSC extraction. Already extracted in ' csc_dir]);
    disp('======================');
    return;
end


%% Extract data for each Tetrode

for ii_Tetrode=(p.TT)
    filename_out_save = [csc_dir, 'CSC_extracted_animal', num2str(p.animal), '_Day', num2str(p.day), '_TT', num2str(ii_Tetrode)];
    
    %% Extract data from the behavioral session
    num_windows = floor(diff(times_microsec_behav_session_2)/Window_size_saved)-1; % Note that we loss ~ 0-1 min
    % of data on this, but I do this because we are taking a window for the
    % filtering twice the size of the window we wish to save and hence if we
    % were to use the maximal number of windown (1 more then what we currently
    % use) we would be using data outside this session.
    %num_windows = floor(diff(times_microsec_behav_session_2)/Window_size_saved);
    
    Num_of_behav_analyzed_files = num_windows;
    
    parfor ii_window = 1:num_windows
        filename_out{ii_window}= [filename_out_save,'_behav_min_num_',num2str(ii_window),'.mat'];
        NLG_filter_CSC_window(filename_out{ii_window},...
            times_microsec_behav_session_2,...
            ii_Tetrode,...
            ii_window,...
            Recording_directory_general,...
            num_windows)
    end
    
end % end looping over Tetrodes

Num_of_analyzed_files.Num_of_behav_analyzed_files = Num_of_behav_analyzed_files;

save ([csc_dir 'Num_of_analyzed_files'], '-struct' ,'Num_of_analyzed_files');