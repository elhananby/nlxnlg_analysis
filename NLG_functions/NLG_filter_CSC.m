function NLG_filter_CSC (p)
global useGPU;
dbstop if error

%% PARAMETERS:
Filter__N_order = 300 ; %  Same as we use for filtering ripples
Filter_FreqsHz = [600 6000]; % Filter for spikes
CSC_Sampling_Rate_Hz = p.nlg.fs;
Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');

%% Define filter
Wn = Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
win_filter = fir1( Filter__N_order, Wn, 'bandpass' ); % Filter parameters are defined above

%% Extract data for each Tetrode
active_channels = find(p.active_channels == 1);

for ii_channel = 1:length(active_channels)
    tic
    fprintf('Fitering Channel %i', ii_channel);
    
    filename_out = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered', active_channels(ii_channel)));
    
    if exist([filename_out '.ncs'], 'file')
         fprintf('\t%.2f\n', toc);
        continue;
    end
    
    filename_CSC_EEG_in = fullfile(Recording_directory_general, ['CSC', num2str(active_channels(ii_channel)-1), '.ncs']);
    
    [Timestamps, Samples, Header] = Nlx2MatCSC( filename_CSC_EEG_in, [1 0 0 0 1], 1, 1, []);

    Samples_to_filter = Samples(:) - mean(Samples(:));
    
    Samples_filtered = FiltFiltM(win_filter, 1, Samples_to_filter);
    
    block_size = 512;
    Samples_save = reshape(Samples_filtered, block_size, length(Samples_filtered)/block_size);
    
    % save ncs
    AppendToFileFlag = 0;
    ExportMode = 1;
    ExportModeVector = [];
    FieldSelectionFlags = [1 0 0 0 1 1];
    
    Mat2NlxCSC([filename_out '.ncs'], AppendToFileFlag, ExportMode, ExportModeVector,...
        FieldSelectionFlags, Timestamps, Samples_save, Header);
    
    fprintf('\t%.2f\n', toc);
end % channels

end