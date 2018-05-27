function Nlg_Flight_detect_spikes_CSC_parfor(p)
dbstop if error
filename_out_temp = fullfile(p.path_dataout, p.datadir_out, 'TT1_cleaned.ntt');
if exist(filename_out_temp, 'file')
    return;
end

%% check if already done spike cleaning

coincidence_window = 500 ; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes

library_file = 'library_of_acceptable_spike_shapes.mat'; % File containing the Library Of Acceptable Spike Shapes

Timestamps_accepted_spikes_TT = [];% Initialize
spikes_TT = [];% Initialize

tetrodes = p.use_tetrodes;
TT_counter = 0;

% Loop over tetrodes, a TEMP file is saved for each tetrode and loaded in
% the next step:
fprintf('Spike Extraction:\n');
for ii_Tetrode = tetrodes 
    fprintf('Processing Tetrode %i...\n', ii_Tetrode);
    filename_out{ii_Tetrode} = Nlg_Flight_detect_spikes_CSC_per_TT_maya(p, ii_Tetrode);  
end %TTs

for ii_Tetrode = tetrodes
    load([filename_out{ii_Tetrode} '.mat']);
    Timestamps_accepted_spikes_TT{ii_Tetrode} = Timestamps_accepted_spikes_per_TT;
    spikes_TT{ii_Tetrode} = spikes_per_TT;
%     delete(filename_out{ii_Tetrode});
    clear Timestamps_accepted_spikes_per_TT spikes_per_TT
end

%%  Read Timestamp data and use Coincidence-Detection across Tetrodes to eliminate artifacts (just as we do in the wired case): --------

idx_coincidence_vec{tetrodes(end)} = [];  % Initialize this variable (for later)

% Find coincidence-detection events = coincidence-detection on a millisecond-scale:
for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{tetrodes(1)}) % Loop over the spikes of the FIRST tetrode
    
    t_spike = Timestamps_accepted_spikes_TT{tetrodes(1)}(ii_spikes);
    idx{tetrodes(1)} = ii_spikes; % Index of this spike
    
    for ii_file = 2:1:length(tetrodes) % Loop over the other tetrodes, finding the coincidnce-detection
        idx{tetrodes(ii_file)} = find( abs( Timestamps_accepted_spikes_TT{tetrodes(ii_file)} - t_spike ) <= coincidence_window ); % THE COINCIDENCE DETECTION
    end
    
    % Check that the coincidence-detection occurred on ALL tetrodes:
    test_variable = 1;
    for ii_file =  tetrodes % Loop over all tetrode files
        if ( isempty( idx{ii_file} ) ) % If there is NO coincidence-detection
            test_variable = 0 ;
        end
    end
    
    if length(tetrodes)>=3 %If there are at least 3 tetrodes
        % Check that the temporal separation between the OTHER two tetrodes also meets the coincidence-detection criterion
        % (this is needed since the spikes on the OTHER two tetrodes may be up to twice-the-time apart, in principle!!! ):
        if ( test_variable == 1 ) % If we passed the first test
            if ( abs( Timestamps_accepted_spikes_TT{tetrodes(2)}(idx{tetrodes(2)}(1)) - Timestamps_accepted_spikes_TT{tetrodes(3)}(idx{tetrodes(3)}(1)) ) > coincidence_window ),
                test_variable = 0 ; % Reset test_variable if the time between spikes on the OTHER two tetrodes is too large
            end
        end
    end
    
    % Save the indexes of the coincidence-detection "spikes" to be REMOVED:
    if ( test_variable == 1 ) % If all indexes exist = there IS a coincidence detection on ALL tetrodes
        for ii_file = tetrodes % Loop over all tetrode files
            idx_coincidence_vec{ii_file} = [ idx_coincidence_vec{ii_file}  idx{ii_file} ];
        end
    end
%     if mod(ii_spikes,1000) == 0
%         disp(['Processing coindicedence detection across TT- ' num2str((ii_spikes/length(Timestamps_accepted_spikes_TT{tetrodes(1)}))*100)])
%     end
%     
end %end looping over spikes for coincidence detection


%% Save the data in Neuralynx NTT files 
for ii_Tetrode = p.use_tetrodes

    spikes = spikes_TT{ii_Tetrode};
    
    %removing spikes because of coincidence detection
    Timestamps_accepted_spikes = Timestamps_accepted_spikes_TT{ii_Tetrode};
    if ~isempty(idx_coincidence_vec)
        spikes = spikes_TT{ii_Tetrode};
        
        %write the spikes that removed due to coincidence detection
        %seperatly:
        Timestamps_accepted_spikes_coincidence = Timestamps_accepted_spikes(idx_coincidence_vec{ii_Tetrode});
        spikes_of_coincidence = spikes(:,:,idx_coincidence_vec{ii_Tetrode});
        
        %remove the coincidence:
        Timestamps_accepted_spikes(idx_coincidence_vec{ii_Tetrode}) = [];
        spikes(:,:,idx_coincidence_vec{ii_Tetrode}) = [];
    end

    FieldSelection = [1 0 1 0 1 1];
    % Timestamps, ~, Cell Numbers, ~, ~, Samples, Header
    CellNumbers = zeros(1, length(Timestamps_accepted_spikes));
    headerfile = 'D:\Scripts\nlx_analysis\NLG_functions\NTT_Header.mat';
    load(headerfile);
    
    if exist([filename_out{ii_Tetrode},'.NTT'], 'file')
        continue;
    end
    
    Mat2NlxSpike([filename_out{ii_Tetrode},'.NTT'], 0, 1, 1, FieldSelection,...
        Timestamps_accepted_spikes, CellNumbers, spikes, Header);
    
end % end looping over tetrodes