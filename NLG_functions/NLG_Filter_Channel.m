function  NLG_Filter_Channel(p)
dbstop if error

%% PARAMETERS:
Filter_N_order = 300 ; %  Same as we use for filtering ripples
Filter_FreqsHz = [600 6000]; % Filter for spikes
Recording_directory_general = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');
CSC_Sampling_Rate_Hz = struct2array(load(fullfile(Recording_directory_general, 'params.mat'), 'fs'));

%% Define filter
Wn = Filter_FreqsHz / ( CSC_Sampling_Rate_Hz / 2); % Passband, Hz (Normalized Wn to half the sampling rate)
win_filter = fir1( Filter_N_order, Wn, 'bandpass' ); % Filter parameters are defined above

%% Extract data for each Tetrode
active_channels = find(p.active_channels == 1);
for ii_channel = 1:length(active_channels)
    tic
    fprintf('Extracting CSC of Day %i from Channel %i', p.day, active_channels(ii_channel));
    
    filename_out = fullfile(Recording_directory_general,...
        sprintf('CSC%i_filtered.dat', active_channels(ii_channel)));
    
    if exist(filename_out, 'file')
        fprintf('\tAlready processed, skipping\n');
        continue;
    end
    
    filename_CSC_EEG_in = fullfile(Recording_directory_general, ['CSC', num2str(active_channels(ii_channel)), '.dat']);
    
    fidIn = fopen(filename_CSC_EEG_in, 'r');
    Samples = fread(fidIn, 'int16');
    fclose(fidIn);
    
    Samples_to_filter = Samples - mean(Samples);
    
    Samples_filtered = FiltFiltM(win_filter, 1, Samples_to_filter);

    fidOut = fopen(filename_out, 'w');
    
    fwrite(fidOut, Samples_filtered, 'int16');
    
    fclose(fidOut);
    
    fprintf('\t%.2f Seconds\n', toc);
end % channels

end