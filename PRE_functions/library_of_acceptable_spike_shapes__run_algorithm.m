function  [ vector_of_sccepted_spikes  vector_of_max_r_values ] = ...
    library_of_acceptable_spike_shapes__run_algorithm(ii_file, library_file, Samples, r_threshold )

% LIBRARY_OF_ACCEPTABLE_SPIKES__RUN_ALGORITHM
%
% SYNTAX:
% [ vector_of_sccepted_spikes  vector_of_max_r_values ] = ...
%      library_of_acceptable_spike_shapes__run_algorithm( library_file, Samples, r_threshold )
%
% PURPOSE:
% This function runs the algorithm for cleaning spike shapes extracted from an Ntt file, based on 
% correlations of the spike's shape with a library of acceptable spikes. For each spike, the
% correlation is computed for the tetrode channel (one out of 4) which has the largest height.
% This correlation is computed for 3 lags (-1, 0, +1), to allow for a shift of +/-1 in the
% position of the waveform's peak -- and the largest of these 3 correlations is taken, and compared 
% to the varabile 'r_threshold', to determine whether to accept or reject this spike.
%
% INPUT PARAMETERS:
%    library_file = mat-file containing a matrix with the acceptable file shapes.
%    Samples = Spike waveforms matrix extracted from an Ntt
%        (tetrode) file using Nlx2MatSpike(); this is the matrix to be cleaned.
%    r_threshold = The algorithm rejects a waveform if it has with correlation below this threshold 
%        with ALL the library's shapes.
%
% OUTPUT PARAMETERS:
%    vector_of_sccepted_spikes = contains '1' for an accepted spike and '0' for a rejected spike.
%    vector_of_max_r_values = maximal correlation (r) value for all the spikes (maximum over all
%        the spikes in the acceptable spikes library).


load('D:\Scripts\nlx_analysis\library_of_acceptable_spike_shapes.mat'); % Load the library of acceptable spike shapes; 
                               % the varaible there should have the name 'library_of_acceptable_spike_shapes'.

vector_of_sccepted_spikes = zeros( 1, size(Samples,3) ) + NaN ; % Initialize
vector_of_max_r_values = zeros( 1, size(Samples,3) ) + NaN ;

% initialize progress bar
m = size(Samples, 3)/100;
fprintf('\nTT%d progress:', ii_file);
fprintf(['\n' repmat('.', 1, 100) '\n\n']);

pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

parfor ii_spike = 1:size(Samples,3) % Loop over spikes extracted from the Ntt file
    
    % progress bar printing
    if mod(ii_spike, floor(m)) == 0
        fprintf('\b|\n');
    end
    
    % Extract the current spike shape on all 4 channels:
    spike_shape_4channels = squeeze( Samples(:,:,ii_spike) ); 
      
    % Choose the channel # for which the spike has the largest height:
    [ stam  idx_channel_max_height ] = max( max( spike_shape_4channels ) ); 
    spike_shape = spike_shape_4channels( :, idx_channel_max_height )' ;
    
    if ( std( spike_shape(2:end-1) ) == 0 ) % If this is a completely FLAT "spike", I cannot compute CORRCOEF, so I will set r = 0
        vector_of_max_r_values( ii_spike ) = 0 ;  % Set r = 0 in this case
        
    else % If this spike DOES have some shape (this is the case basically for ALL the recorded waveforms)
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag 0 :
          xxx_lags = 2 : 31 ; % Use CENTRAL 30 points
          % run for the current spike shape in comparison to all acceptable
          % spike shpes.
          ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
          rrr_vec_lag_0 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag (+1) :
          xxx_lags = 1 : 30 ; % Use FIRST 30 points (RIGHT shift of the "Acceptable Spike Shapes matrix")
          ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
          rrr_vec_lag_plus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Compute the correlation coefficients with all the acceptable spike shapes -- lag (-1) :
          xxx_lags = 3 : 32 ; % Use LAST 30 points (LEFT shift of the "Acceptable Spike Shapes matrix")
          ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
          rrr_vec_lag_minus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
        % Save the MAXIMAL r value -- use the maximal correlation among the 3 lags (-1,0,+1):
          vector_of_max_r_values( ii_spike ) = max( [ rrr_vec_lag_0  rrr_vec_lag_plus1  rrr_vec_lag_minus1 ] );
    end
            
    % Determine if this spike should be Accepted (set value to '1') or Rejected (set value to '0'):
    vector_of_sccepted_spikes( ii_spike ) = ( vector_of_max_r_values( ii_spike )  >=  r_threshold ); 
        % Accept the spike shape ('1') if its correlation with ANY of the acceptable shapes 
        % is >= r_threshold ; else, reject the spike ('0').
    
end % End "Loop over spikes extracted from the Ntt file"

% --- End ---
