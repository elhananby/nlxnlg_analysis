function NLG_PRE_create_inclusion_list_of_neurons(P,C)


% For each Day-Tetrode combination, I plot a figure showing:
% (1) Activity pattern of all the neurons recorded on this tetrode.
% (2) Spike Waveforms for the neurons (for all 4 TT channels).
% (3) ISI Histograms for the neurons.
% (4) 2-D and 3-D cluster plots (energy plots for all the spikes) for the data from this tetrode, colored by-neuron.
% *** This script also saves these figures into D:\Michael\Data\Expdata_Processed\Combined_data\ .
%
% The output of this script is the mat-file:  D:\Michael\Data\Expdata_Processed\Combined_data\data_inclusion_list_of_neurons.mat ,
% which contains:
% (1) "file_list_allunits" :  The filenames of ALL the neurons (with full
% path).
% (2) "file_list_associated_VT_files_allunits" :  VT Filename (actually a mat-file) associated with these Spike files.
% (3) "spike_width_ms_allunits" :  Spike width, peak-to-trough (ms).
% (4) "average_firing_rate_Hz_allunits" :  Mean firing rate, averaged over ALL 5 (or 3) sessions.
% (5) "number_of_spikes_in_every_session_allunits" :  Number of spikes in every session (NaN if the session was not run).
% (6) "inclusion_list_criteria_comments" :  The Criteria used for making the assignments below.
% (7) "is_unit_is_pyramidal_is_active_allunits" :  [1 1 1] means a cell that is: [acceptable single unit, pyramidal cell, behaviorally-active].
% FOR THE LAST VARIABLE, I WRITE THESE '1' AND '0' VALUES ***MANUALLY*** IN THE BATCH-FILE.
%
% As a byproduct, this program also saves, for each Day-Tetrode combination, a COMBINED Ntt file,
% containing the spikes from ALL my spike-sorted neurons, color-coded by-neuron (this is useful for
% doing quick visualization of the data using SpikeSort3D.exe). This Ntt file is saved in the
% same directory in D:\Michael\Data\Expdata_Processed\ , where my individual spike-sorted Ntt files are stored.
%


%-----------------
% Dori D, adapted from Michael Yartsev, Maya G
%-----------------
%


% ---------  General Parameters: ----------------------------

p = P(1);
p.path_dataout='D:\Ayelet\Data\Data_Nlg_Proc\Homing_vector\'
filename_out = fullfile(p.path_dataout,'inclusion_list_of_neurons.mat') ; % mat-file for saving the inclusion list
dir_figs_out = fullfile(p.path_dataout,'Homing_vector\Inclusion_List_figures'); ; % Directory for saving the figures (1 figure for each Day-Tetrode combination)
amplitude_scale_bar_val = 100; %The value of the amplitude bar we want to be displayed in the figures.

inclusion_list_criteria_comments = [ ...
    '(1) ACCEPTNING A NEURON INTO THE DATABASE:                                                                          '; ...
    '    ***MANUAL*** SELECTION, written by me MANUALLY in the batch-file -- based on following CRITERIA:                '; ...
    '    (a) TOTAL Spike Count in ALL my 5 sessions (or 3 sessions) >= 100 spikes.                                       '; ...
    '    (b) Spike shape should be neural-like = HAVE Afterhyperpolarization + NOT have a funky shape.                   '; ...
    '    (c) ISI Histogram should have no intervals < 2 ms  (i.e. this is NOT a mixture of 2 cells).                     '; ...
    '    (d) The SC Threshold should be clearly BELOW THE SMALLEST WAVEFORMS of the unit: i.e. this cluster was not cut. '; ...
    '(2) ASSIGNING A CELL AS A PYRAMIDAL CELL:                                                                           '; ...
    '    ***MANUAL*** SELECTION, written by me MANUALLY in the batch-file -- based on following CRITERIA:                '; ...
    '    (a) A pyramidal Shape = narrow high peak followed by a long shallow trough (afterhyperpolarization).            '; ...
    '    (b) Average Firing Rate over all 5 sessions <= 5 Sp/s  (the number is from: Gothard et al., J. Neurosci. 1996). '; ...
    '    (c) ISI Histogram has a clear double-peak, demonstrating COMPLEX SPIKES (BURSTS).                               '; ...
    '(3) ASSIGNING A CELL AS A "SILENT CELL" (DURING BEHAVIOR) VERSUS "ACTIVE CELL":                                     '; ...
    '    ***MANUAL*** SELECTION, written by me MANUALLY in the batch-file -- based on following CRITERIA:                '; ...
    '    (a) Spike Count in AT LEAST ONE of my 2 behavioral sessions >= 40 spikes  (40 spikes / 20 min = 2 spikes/min).  '
    ];

color_list = [ 0 0 1 ; ... % List of 18 colors for plotting spikes of different units
    0 1 0 ; ...
    1 0 0 ; ...
    0.5 0.5 0 ; ...
    0.5 0 0.5 ; ...
    0 0.5 0.5 ; ...
    0.9 0.1 0.5 ; ...
    0.9 0.5 0.1 ; ...
    0.1 0.9 0.5 ; ...
    0.1 0.5 0.9 ; ...
    0.5 0.9 0.1 ; ...
    0.5 0.1 0.9 ; ...
    0.3 0.7 0.7 ; ...
    0.7 0.3 0.7 ; ...
    0.7 0.7 0.3 ; ...
    0.3 0.3 0.7 ; ...
    0.3 0.7 0.3 ; ...
    0.7 0.3 0.3 ] ;

max_neurons_to_display = 16; % maximum number of neurons to display in plot
% (per tetrode per day)
% --------------------------------------------------------------


% eval(['cd ', data_dir]);


       
% ======== Saving the Data for all the neurons: ===========

inclusion_list_index = 0; % index into inclusion list

% enumerate on tetrode-day combination

figure ;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[3 2 29 22]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0.2 0.2 0 0]);

for ii_experiment = 1:length(P)
    %% Ayelet- find behavioral session timestemp
% find the string of start/end session in the events strings:
        find_end_session_1_event=strfind(P(1, ii_experiment).Nlg_EventStrings,'End session 1');
        find_start_session_2_event=strfind(P(1, ii_experiment).Nlg_EventStrings,'Start session 2');
         
         % adding session 3 if exist:
        if ~isempty(strfind(P(1, ii_experiment).Nlg_EventStrings,'Start session 3'))
            find_end_session_2_event=strfind(P(1, ii_experiment).Nlg_EventStrings,'End session 2');
            find_start_session_3_event=strfind(P(1, ii_experiment).Nlg_EventStrings,'Start session 3');
            end_session_2_event=find(~cellfun(@isempty,find_end_session_2_event));
            start_session_3_event=find(~cellfun(@isempty,find_start_session_3_event)); 
            end_session_2_timestamp=P(1, ii_experiment).Nlg_EventTimestamps (end_session_2_event);
        start_session_3_timestamp=P(1, ii_experiment).Nlg_EventTimestamps (start_session_3_event);
        end
        
        end_session_1_event=find(~cellfun(@isempty,find_end_session_1_event));
        start_session_2_event=find(~cellfun(@isempty,find_start_session_2_event));
        %find the event timestamp:
        end_session_1_timestamp=P(1, ii_experiment).Nlg_EventTimestamps (end_session_1_event);
        start_session_2_timestamp=P(1, ii_experiment).Nlg_EventTimestamps (start_session_2_event);
        
    p = P(ii_experiment); % copy experiment data to structure "p"
    bat = p.bat;
    day = p.day;
    data_dir = sprintf('%s\\%s\\%d\\',p.path_dataout,p.path_year_bat,p.day);
    session_names = {p.S.session};
    nsessions = p.nsessions;
    
    
    % enumerate over tetrode
    
    for tetrode = p.TT
        
        Timestamps_combined = []; % Initialize -- here I will save the data to be written to the Ntt file
        ChanNum_combined = [];
        CellNumbersSpikeSorting_combined = [];
        NumValidSamples_combined = [];
        Samples_combined = [];
        
        % find nominal depth of current tetrode
        
        tetrode_index = find(ismember(p.TT,tetrode));
        nominal_depth = p.depth(tetrode_index);
        
        % enumerate over cells within that day/tetrode combination
        
        cells_index = find(ismember([C.bat],bat) & ...
            ismember([C.day],day) & ...
            ismember([C.TT],tetrode));
        
        if isempty(cells_index)
            continue;
        end
        
        
        first_index_for_day = inclusion_list_index+1;
        
        for cell_index = cells_index
            
            inclusion_list_index = inclusion_list_index + 1; % index into inclusion cell list
            
            c = C(cell_index);
            
            cell_id = c.cell_id;
            
            
            % --- (1) "file_list_allunits" : ---
            
            
            file_list_allunits{ inclusion_list_index } = ...
                sprintf('%s\\bat%d_Day%d_TT%d_SS_%02d.ntt', ...
                data_dir,p.bat,p.day,tetrode,cell_id);
            
            file_list_in_all_threshold_crossing_un_sorted_all_units{ inclusion_list_index } = ...
                sprintf('%s\\bat%d_Day%d_TT%d.ntt', ...
                data_dir,p.bat,p.day,tetrode); % This will be used for the multi-unit autocorrelogram in later processing
            
            file_names_in{ inclusion_list_index } = ...
                sprintf('bat%d_Day%d_TT%d_SS_%02d.ntt', ...
                p.bat,p.day,tetrode,cell_id);
            
            experiment_parameters_allunits{ inclusion_list_index } = p;
            cell_parameters_allunits{ inclusion_list_index } = c;
            
            % also keep a file for every session seperately
            for nses = 1:nsessions
                session = p.S(nses).session;
                file_list_session_allunits{ inclusion_list_index }{nses} = ...
                    sprintf('%s\\bat%d_Day%d_%d_%s_TT%d_SS_%02d.ntt', ...
                    data_dir,p.bat,p.day,nses,session,tetrode,cell_id);
                
                file_list_session_in_all_threshold_crossing_un_sorted_all_units{ inclusion_list_index }{nses} = ...
                    sprintf('%s\\bat%d_Day%d_%d_%s_TT%d.ntt', ...
                    data_dir,p.bat,p.day,nses,session,tetrode);
                
                session_name_list_all{ inclusion_list_index }{nses} = session;
                
            end
            
            % --- Collect data for the "combined" Ntt file: ---
            
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( file_list_allunits{ inclusion_list_index }, ...
                [1 1 1 1 1], 1, 1, [] ) ; % Read all the data
            
            Timestamps_combined = [Timestamps_combined  Timestamps];
            ChanNum_combined = [ChanNum_combined  ChanNum];
            CellNumbersSpikeSorting_combined = [CellNumbersSpikeSorting_combined  CellNumbersSpikeSorting];
            NumValidSamples_combined = cat( 2, NumValidSamples_combined, NumValidSamples); % Concatenating 2-D arrays
            Samples_combined = cat( 3, Samples_combined, Samples ); % Concatenating 3-D arrays
            
            % --- (2) "file_list_associated_VT_files_allunits" :  VT Filename (actually a mat-file) associated with these Spike files. ---
            
            for nses = 1:nsessions
                %                 VT_filename{nses} = ['VT_extracted_bat',num2str(p.bat),'_Day',num2str(p.day) ...
                %                     '_' num2str(nses) '_' session_names{nses}];
                VT_filename{nses} = ['VT_flight_extracted_bat',num2str(p.bat)];
                file_list_associated_VT_files_allunits{inclusion_list_index}{nses} = ...
                    fullfile(data_dir,VT_filename{nses});% I will save this filename as well
            end
            
            % --- (3) "spike_width_ms_allunits" (peak-to-trough width) -- computed over SPLINED ("UP-SAMPLED") waveform: ---
            
            keyword = '-SamplingFrequency '; % notice the trailing blank
            header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)')); % usually 13
            %             SampleFrequency_Hz = str2num( NlxHeader{header_index}(length(keyword)+1:end) ); % Sample frequency on SC channels
            [stam  idx_maximal_channel] = max( max( squeeze( mean( Samples, 3 ) ) ) ); % TT Channel with the largest average-spike
            waveform = mean( squeeze( Samples( :, idx_maximal_channel, : ) )' ); % Waveform SHAPE
            interp_Multiply_Sample_Frequency_by_Factor = 10 ; % For spline interpolation, use a sampling rate that is faster by this factor
            xxx_interp = (1+1/interp_Multiply_Sample_Frequency_by_Factor) : (1/interp_Multiply_Sample_Frequency_by_Factor) : 32 ; % x-data-points for spline interpolation of the Waveform
            waveform_interp = spline( 1:32, waveform, xxx_interp ); % Computing spike width over SPLINED ("UP-SAMPLED") waveform
            [stam  idx_peak] = max( waveform_interp );
            [stam  idx_trough] = min( waveform_interp );
            %             spike_width_ms_allunits( inclusion_list_index ) = ... % Save this data
            %                 (idx_trough - idx_peak)/ SampleFrequency_Hz / interp_Multiply_Sample_Frequency_by_Factor * 1000 ;
            %
            
            % --- (4) "average_firing_rate_Hz_allunits" (averaged over ALL sessions): ---
            
            start_times = [p.S.start_time]; %it was *1000 but I think it was relevant for an old code..
            end_times = [p.S.end_time];%it was *1000 but I think it was relevant for an old code..
            total_duration_of_all_sessions = sum(end_times-start_times);
            average_firing_rate_Hz_allunits( inclusion_list_index ) = ... % Save this data
                length( Timestamps ) / total_duration_of_all_sessions * 10^6 ;
            
            
            % --- (5) "number_of_spikes_in_every_session_allunits" (for ALL sessions): ---
            
            nsessions = length(p.S);
             start_times = [p.S.start_time]; %it was *1000 but I think it was relevant for an old code..
            end_times = [p.S.end_time];%it was *1000 but I think it was relevant for an old code..

            for nses = 1:nsessions
                number_of_spikes_in_every_session_allunits( inclusion_list_index, nses ) = ...
                    sum( Timestamps >= start_times(nses) & Timestamps < end_times(nses) );
                session_start_end_time{inclusion_list_index,nses}=[start_times(nses),end_times(nses) ]; 
            end
            
            % --- (6) "inclusion_list_criteria_comments": ---
            
            inclusion_list_criteria_comments = inclusion_list_criteria_comments ; % This is the same set of commetns for ALL the cells
            
            
            % --- (7) "is_unit_is_pyramidal_is_active_allunits" : ---
            
            % [1 1 1] means a cell that is: [acceptable single unit, pyramidal cell, behaviorally-active].
            
            is_unit_is_pyramidal_is_active_allunits(inclusion_list_index ,:) = ...
                [c.single_unit c.pyramidal c.behaviorally_active ];
            
            % --- (8) "Peak-to-Peak amplitude of the average waveform" : ---
            
            keyword = '-ADBitVolts '; % notice the trailing blank
            header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)')); % usually 15
            NlxHeader_ADBitVolts = str2num( NlxHeader{header_index}(length(keyword)+1:end) ); % Sample frequency on SC channels
            NlxHeader_ADBitVolts = NlxHeader_ADBitVolts(1);
            peak_to_peak_amp_micro_volts(inclusion_list_index)=0; % initialize
            for ii_channel = 1:4,
                spike_shapes = squeeze( Samples(:,ii_channel,:) );
                % Compute the peak-to-peak value of the average spike wavewform.
                peak_to_peak_amp_micro_volts_temp = ( max(mean(spike_shapes,2))-min(mean(spike_shapes,2)))*10^6*NlxHeader_ADBitVolts;
                if (peak_to_peak_amp_micro_volts_temp>peak_to_peak_amp_micro_volts(inclusion_list_index)) % Find the peak to peak on the maximal channel
                    peak_to_peak_amp_micro_volts(inclusion_list_index) = peak_to_peak_amp_micro_volts_temp;
                else end
            end
            
            % --- (9) "Times of Beheav Session" : ---
            
            nsessions = length(p.S);
             start_times = [p.S.start_time]; %it was *1000 but I think it was relevant for an old code..
            end_times = [p.S.end_time];%it was *1000 but I think it was relevant for an old code..

            behavior_sessions = p.use_for_analysis;
            behav_times = [start_times(behavior_sessions);end_times(behavior_sessions)]';
            time_of_behav_session_all_units{ inclusion_list_index } = behav_times;
            
            
            %             % --- (10) "Voltage Threshold" : ---
            %
            %             keyword = '-ThreshVal '; % notice the trailing blank
            %             header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)'));
            %             thresh_vals = str2num( NlxHeader{header_index}(length(keyword)+1:end));
            %             NlxHeader_Threshold(inclusion_list_index ) = max(thresh_vals); % A/D bits-to-volts conversion
            
            % --- (11) "nominal depth" : ---
            
            nominal_depth_all_units( inclusion_list_index ) = ...
                nominal_depth;
            
            Sync_by_LED( inclusion_list_index ) = p.Sync_by_LED;
            
        end % for ii_cell ... (loop over cells for this bat-day-tetrode combination)
        
        
        % ========= PLOT A SEPARATE FIGURE FOR EVERY BAT-DAY-TETRODE: ===========
        
        clf; %clean figure
        
        
        % -------- Plot Activity Patterns & Mark the Sessions: --------
        
        axes( 'position', [0.04 0.06 0.25 0.90]) ; hold on;
        
        % Set some graphical parameters:
        xticks = sort([start_times end_times]); % Put xticks at these times
        xticks_labels = round( xticks / 10^6); % Put ROUNDED-NUMBERS for labels, and put them in sec
        max_neurons_to_display = 16;
        ylimits = [-1 max_neurons_to_display+1]; % This well-covers the maximal number of neurons in any of my recording days-tetrodes
        xlimits(1) = min( start_times );
        xlimits(2) = max( end_times );
        xlimits = [xlimits(1)-300*10^6  xlimits(2)+300*10^6]; % Give some ore space on both sides
        
        number_of_cells_today = inclusion_list_index - first_index_for_day + 1;
        
        
        idx_units_on_this_tetrode = [ first_index_for_day:inclusion_list_index ]; % The indexes of the cells on this tetrode, from the "allunits" list
        file_list_in = file_list_allunits( idx_units_on_this_tetrode );
        
        for ii_cell = 1:number_of_cells_today % Loop over the Cells (for this bat-day-tetrode combination)
            
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( file_list_in{ ii_cell }, [1 1 1 1 1], 1, 1, [] ) ; % Read all the data
            
            % Determine the color-of-dots to plot for the activity pattern
            % (CYAN = [0 1 0] = non-acceptable cell; BLUE = [1 1 0] = Acceptable, but Non-Active ; BLACK = [1 1 1] = Active cell):
            
            color_dots = 'c' ; % Default (for "is_unit" = '0' = non-acceptable cell)
            if ( is_unit_is_pyramidal_is_active_allunits( idx_units_on_this_tetrode( ii_cell ), 1 ) == 1 ), % If "is_unit" = '1' = acceptable cell
                color_dots = 'y' ; % yellow
            end
            if ( is_unit_is_pyramidal_is_active_allunits( idx_units_on_this_tetrode( ii_cell ), 3 ) == 1 ), % If "is_active" = '1' = Active cell
                color_dots = 'g' ; % green
            end
            
            % Plot the Activity Pattern:
            
            plot( Timestamps, ii_cell, '.', 'markersize', 5, 'color', color_dots );hold on;
            plot(start_session_2_timestamp,ii_cell,'*','markersize', 5, 'color','k');
           
            % Mark the Sleep Sessions (Gray lines) and Behavioral Sessions (Red lines):
            
            for nses = 1:nsessions
                
                if ismember(nses,p.use_for_analysis) % behavioral session
                    color_spec = 'r';
                else % sleep sessopm
                    color_spec = [0.7 0.7 0.7];
                end
                hhline = line( [start_times(nses) end_times(nses)], [-0.5 -0.5], 'color', color_spec, 'linewidth', 8);
                nspikes = sum( Timestamps >= start_times(nses) & Timestamps < end_times(nses) );
                text( start_times(nses), ii_cell + 0.3, ['n=', num2str(nspikes)], 'fontsize', 7 );
                text( start_times(nses), -0.8, session_names{nses}, 'fontsize', 7);
            end % nses
            
        end % End Loop over the cells (for this bat-day-tetrdoe combination)
        
        % Set, etc:
        xlabel('Time (s)');
        ylabel('Neuron #');
        title( data_dir, 'interpreter', 'none', 'fontweight', 'bold' );
        set( gca, 'ylim', ylimits, 'ytick', [1:length(file_list_in)], 'xlim', xlimits, ...
            'xtick', xticks, 'xticklabel', xticks_labels, 'tickdir', 'out', 'ticklength', [0.003 0.003] );
        
        
        
        % -------- Plot Waveforms: --------
        
        % Set some graphical parameters:
        horizontal_separation_waveforms = 45 ; % Horizontal separation for plotting spike waveforms (separation between the TT's 4 channels)
        xlimits = [0   40+horizontal_separation_waveforms*(4-1)]; % The right edge is the right-end of the waveform in ch. 4
        
        for ii_cell = 1:number_of_cells_today % Loop over the Cells :
            
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( file_list_in{ ii_cell }, [1 1 1 1 1], 1, 1, [] ) ; % Read all the data
            
            ylimits_this_unit = [ min(Samples(:))  max(Samples(:)) ];
            
            axes( 'position', [0.30  0.135+(ii_cell-1)*0.0495  0.25  0.0595] ); hold on;
            
            % --- Plotting the Spike Shapes ---
            keyword = '-ADBitVolts '; % notice the trailing blank
            header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)')); % usually 15
            NlxHeader_ADBitVolts = str2num( NlxHeader{header_index}(length(keyword)+1:end) ); % Sample frequency on SC channels
            NlxHeader_ADBitVolts = NlxHeader_ADBitVolts(1); % A/D bits-to-volts conversion, for CHANNEL #1
            
            for ii_channel = 1:4,
                spike_shapes = squeeze( Samples(:,ii_channel,:) );
                plot( (1:32) + horizontal_separation_waveforms*(ii_channel-1), spike_shapes, ...
                    'color', color_list(ii_cell,:), 'linewidth', 0.25 ); % I use the variable 'color_list', defined ABOVE
                % Compute the peak-to-peak value of the average spike wavewform.
                %         peak_to_peak_amp_micro_volts_temp = ( max(mean(spike_shapes,2))-min(mean(spike_shapes,2)))*10^6*NlxHeader_ADBitVolts;
                %         if (peak_to_peak_amp_micro_volts_temp>peak_to_peak_amp_micro_volts) % Find the peak to peak on the maximal channel
                %             peak_to_peak_amp_micro_volts = peak_to_peak_amp_micro_volts_temp;
                %         else end
            end
            
            % --- Write the cell's parameters: ---
            
            
            text( xlimits(2)+diff(xlimits)/20, 0, ...
                {
                %['Width (ms) = ', ...
                %                 num2str( spike_width_ms_allunits( idx_units_on_this_tetrode( ii_cell ) ), 2 )], ...
                ['Rate (sp/s) = ', ...
                num2str( average_firing_rate_Hz_allunits( idx_units_on_this_tetrode( ii_cell ) ), 2 )], ...
                ['p-p','[','\muv',']',' = ' ,num2str( peak_to_peak_amp_micro_volts( idx_units_on_this_tetrode( ii_cell ) ) )], }, ...
                'fontsize', 7 );
            %Note that the peak to peak value is on the maximal channel!
            
            
            % [1 1 1] means a cell that is: [acceptable single unit, pyramidal cell, behaviorally-active].
            
            %             keyword = '-ThreshVal '; % notice the trailing blank
            %             header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)'));
            %             thresh_vals = str2num( NlxHeader{header_index}(length(keyword)+1:end));
            %             thresh_val = max(thresh_vals);
            %
            if ( ii_cell == length(file_list_in) ), % For the last cell only,
                text( horizontal_separation_waveforms, ylimits_this_unit(2)+diff(ylimits_this_unit)*1.1, ...
                    {
                    %['Max. Threshold on ALL channels = ', num2str(thresh_val) ,'\muv'],...
                    ['Amplitude scale bar = ', num2str(amplitude_scale_bar_val),'\muv'],...
                    ['CYAN dots (at left) = non-acceptable single unit'], ...
                    ['BLUE = acceptable but non-active cell'], ...
                    ['BLACK = Active cell'], ...
                    ['[1 1 1] = [Single unit, Pyramidal, Behaviorally-active]']}, ...
                    'fontsize', 7, 'fontweight', 'bold');
            end
            
            % Marking the SC Threshold (which , in my experiments, IS IDENTICAL FOR ALL TT CHANNELS!!):
            
            %             NlxHeader_ThreshVal= thresh_val; % Threshold, for CHANNEL #1 (microvolts)
            
            keyword = '-ADBitVolts '; % notice the trailing blank
            header_index = find(~cellfun('isempty',strfind(NlxHeader,keyword)')); % usually 15
            NlxHeader_ADBitVolts = str2num( NlxHeader{header_index}(length(keyword)+1:end) ); % Sample frequency on SC channels
            NlxHeader_ADBitVolts = NlxHeader_ADBitVolts(1); % A/D bits-to-volts conversion, for CHANNEL #1
            %
            %             SC_threshold_in_samples = NlxHeader_ThreshVal / NlxHeader_ADBitVolts / 10^6 ; % Now I can plot the SC Threshold together with the Waveforms
            %             line( [xlimits(1) xlimits(2)-5], [1 1]*SC_threshold_in_samples, 'color', 'red' );
            
            % Marking an amplitude bar which will be the SC Threshold (which , in my experiments, IS IDENTICAL
            % FOR ALL TT CHANNELS!!): We will creat an amplitude bar of 100
            % microvolts.
            
            line( [xlimits(2)-2,xlimits(2)-2], [0 amplitude_scale_bar_val / NlxHeader_ADBitVolts / 10^6], 'color', 'black','linewidth', 2 );
            
            % Set:
            
            set( gca, 'ylim', ylimits_this_unit, 'xlim', xlimits, 'visible', 'off' );
            
        end % End Loop over the cells
        
        
        
        
        % -------- Plot ISI Histograms: --------
        
        % Set some graphical parameters:
        xlimits = [0 5]; % from 1 ms (log10(1)=0) to 100 sec (log10(10^5)=5)
        nbins = 50; % Number of bins in the histogram
        xticklabels = [2 10 1000 100000]; % xticks (ms)
        xticks = log10( xticklabels ); % log10( xticks (ms) )
        
        for ii_cell = 1:number_of_cells_today % Loop over the Cells :
            
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( file_list_in{ ii_cell }, [1 1 1 1 1], 1, 1, [] ) ; % Read all the data
            
            axes( 'position', [0.67  0.135+(ii_cell-1)*0.0495  0.08  0.046] ); hold on;
            
            ISI_ms = diff( Timestamps ) / 10^3 ; % ms
            log10_ISI_ms = log10( ISI_ms );
            
            hist( log10_ISI_ms, nbins );
            
            set( gca, 'xlim', xlimits, 'xtick', xticks, 'xticklabel',[], 'ytick', [], ...
                'ticklength', [0.045 0.045], 'tickdir', 'out' );
            if (ii_cell==1), % For the first cell only
                set( gca, 'xticklabel', xticklabels );
                xlabel('ISI (ms)');
            end
            
        end % End Loop over the cells
        
        
        
        
        % ----- Plot Cluster Displays (Spike-Energy Displays), in 2-D and 3-D: -----
        
        % Initialize the axes:
        hh_1_2 = axes( 'position', [0.82  (0.02 + 0 * 0.12)  0.10  0.11] ); hold on; % Axes for 2-D Plots
        hh_1_3 = axes( 'position', [0.82  (0.02 + 1 * 0.12)  0.10  0.11] ); hold on;
        hh_1_4 = axes( 'position', [0.82  (0.02 + 2 * 0.12)  0.10  0.11] ); hold on;
        hh_2_3 = axes( 'position', [0.82  (0.02 + 3 * 0.12)  0.10  0.11] ); hold on;
        hh_2_4 = axes( 'position', [0.82  (0.02 + 4 * 0.12)  0.10  0.11] ); hold on;
        hh_3_4 = axes( 'position', [0.82  (0.02 + 5 * 0.12)  0.10  0.11] ); hold on;
        hhh_1_2_3 = axes( 'position', [0.82  (0.02 + 6 * 0.12)+0.02  0.10  0.11] ); hold on; % Axes for 3-D Plots
        hhh_2_3_4 = axes( 'position', [0.82  (0.02 + 7 * 0.12)+0.02  0.10  0.11] ); hold on;
        
        
        
        max_energy = 0 ; % Initialize
        
        for ii_cell = 1:number_of_cells_today, % Loop over the Cells :
            
            [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
                Nlx2MatSpike( file_list_in{ ii_cell }, [1 1 1 1 1], 1, 1, [] ) ; % Read all the data
            
            spike_energy_all4channels = []; % Initialize
            for ii_channel = 1:4,
                spike_shapes = squeeze( Samples(:,ii_channel,:) );
                spike_energy_all4channels(:,ii_channel) = sqrt( sum( spike_shapes.^2 ) )';
            end
            
            % --- Plot 2-D Energy-Plots (2-D clusters): ---
            
            axes( hh_1_2 );
            hold on
            xdata = spike_energy_all4channels(:,1) ;
            ydata = spike_energy_all4channels(:,2) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 1', 'Y: Energy 2'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            axes( hh_1_3 );
            hold on
            xdata = spike_energy_all4channels(:,1) ;
            ydata = spike_energy_all4channels(:,3) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 1', 'Y: Energy 3'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            axes( hh_1_4 );
            hold on
            xdata = spike_energy_all4channels(:,1) ;
            ydata = spike_energy_all4channels(:,4) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 1', 'Y: Energy 4'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            axes( hh_2_3 );
            hold on
            xdata = spike_energy_all4channels(:,2) ;
            ydata = spike_energy_all4channels(:,3) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 2', 'Y: Energy 3'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            axes( hh_2_4 );
            hold on
            xdata = spike_energy_all4channels(:,2) ;
            ydata = spike_energy_all4channels(:,4) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 2', 'Y: Energy 4'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            axes( hh_3_4 );
            hold on
            xdata = spike_energy_all4channels(:,3) ;
            ydata = spike_energy_all4channels(:,4) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            plot( xdata, ydata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2), ylimits(2)/2, ...
                    {'X: Energy 3', 'Y: Energy 4'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'xlim', xlimits, 'ylim', ylimits );
            
            % --- Plot 3-D Energy-Plots (3-D clusters): ---
            
            axes( hhh_1_2_3 );
            hold on
            xdata = spike_energy_all4channels(:,1) ;
            ydata = spike_energy_all4channels(:,2) ;
            zdata = spike_energy_all4channels(:,3) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            zlimits = [0  max_energy];
            plot3( xdata, ydata, zdata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            view ( -45, 45 ) ; % Specify the 3-D viewing angle
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2)*1.5, ylimits(2)/2, ...
                    {'X: Energy 1', 'Y: Energy 2', 'Z: Energy 3'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'ztick', [], ...
                'xlim', xlimits, 'ylim', ylimits, 'zlim', zlimits );
            axes( hhh_2_3_4 );
            hold on
            xdata = spike_energy_all4channels(:,2) ;
            ydata = spike_energy_all4channels(:,3) ;
            zdata = spike_energy_all4channels(:,4) ;
            max_energy = max([max_energy ; max(spike_energy_all4channels(:))]); % At the end  of the looping-over-the-cells, this will give the Max over ALL units
            xlimits = [0  max_energy];
            ylimits = [0  max_energy];
            zlimits = [0  max_energy];
            plot3( xdata, ydata, zdata, '.', 'markersize', 2, 'color', color_list(ii_cell,:) ); % I use the variable 'color_list', defined ABOVE
            view ( -45, 45 ) ; % Specify the 3-D viewing angle
            if ( ii_cell == length(file_list_in) ), % Write the text for the Last cell only (when max_energy is finalized)
                text( xlimits(2)*1.5, ylimits(2)/2, ...
                    {'X: Energy 2', 'Y: Energy 3', 'Z: Energy 4'}, ...
                    'fontsize', 7 );
            end
            set( gca, 'xtick', []', 'ytick', [], 'ztick', [], ...
                'xlim', xlimits, 'ylim', ylimits, 'zlim', zlimits );
            
            
        end
        
        
        
        
        
        
        % ======= SAVE FILES: ==============
        
        
        % --- Save the Ntt file with all the spikes for all the units, colored (for quick look in SpikeSort3D.exe): ---
        
        [stam, idx_sort] = sort( Timestamps_combined ); % Sort the data by ascending Timestamps
        Timestamps_combined = Timestamps_combined( idx_sort );
        ChanNum_combined = ChanNum_combined( idx_sort );
        CellNumbersSpikeSorting_combined = CellNumbersSpikeSorting_combined( idx_sort );
        NumValidSamples_combined = NumValidSamples_combined( :, idx_sort );
        Samples_combined = Samples_combined( :, :, idx_sort );
        NlxHeader_combined = NlxHeader ; % This is the same variable for all the units
        
        filname_ntt_out = [ data_dir, '\', file_names_in{1}(1:end-7), '_combined.ntt' ];
        
        Mat2NlxTT( filname_ntt_out, 0, 1, 1, length(Timestamps_combined), [1 1 1 1 1 1], ...
            Timestamps_combined, ChanNum_combined, CellNumbersSpikeSorting_combined, NumValidSamples_combined, Samples_combined, NlxHeader_combined );
        
        
        %--- Save the Figure into a JPG file: ---
        if ~exist(dir_figs_out)
            mkdir(dir_figs_out);
        end;
        filname_fig_out = [ dir_figs_out, '\fig_InclusionList_', file_names_in{inclusion_list_index}(1:end-10), '.jpg' ];
        
        eval( ['print -f', num2str(gcf), ' ', filname_fig_out, ' -djpeg'] );
        
        close all ; % CLOSE THE FIGURE (TOO MANY FIGURES MAY BE DIFFICULT FOR MY COMPUTER):
        
    end % for tetrode ...
    
end % for ii_experiment ...

% add some more parameters to inclusion list from the P and C input structures

bat_vec = []; day_vec = []; TT_vec = []; cell_id_vec = [];
bat_name_list = {}; behavior_sessions_list = {};

ncells = inclusion_list_index;
for inclusion_list_index = 1:ncells
    experiment_parameters = experiment_parameters_allunits{inclusion_list_index};
    cell_parameters = cell_parameters_allunits{inclusion_list_index};
    bat_vec(end+1) = cell_parameters.bat;
    day_vec(end+1) = cell_parameters.day;
    TT_vec(end+1) = cell_parameters.TT;  % tetrode
    cell_id_vec(end+1) = cell_parameters.cell_id; % cell number, according to SpikeSort3D classification
    bat_name_list{end+1} = experiment_parameters.bat_name;
    behavior_sessions_list{end+1} = experiment_parameters.use_for_analysis;
end

% --- Save the MAT-file with all the important variables: ---

inclusion_list.bat = bat_vec;
inclusion_list.day = day_vec;
inclusion_list.TT = TT_vec;
inclusion_list.cell_id = cell_id_vec;
inclusion_list.bat_name = bat_name_list;
inclusion_list.file_list_allunits = file_list_allunits' ; % Save all the important variables in one struct
inclusion_list.time_of_behav_session_all_units = time_of_behav_session_all_units;
inclusion_list.file_list_in_all_threshold_crossing_un_sorted_all_units = file_list_in_all_threshold_crossing_un_sorted_all_units' ; % Save all the important variables in one struct
inclusion_list.file_list_associated_VT_files_allunits = file_list_associated_VT_files_allunits' ;
% inclusion_list.spike_width_ms_allunits = spike_width_ms_allunits' ;
inclusion_list.average_firing_rate_Hz_allunits = average_firing_rate_Hz_allunits' ;
inclusion_list.number_of_spikes_in_every_session_allunits = number_of_spikes_in_every_session_allunits ;
inclusion_list.inclusion_list_criteria_comments = inclusion_list_criteria_comments ;
inclusion_list.is_unit_is_pyramidal_is_active_allunits = is_unit_is_pyramidal_is_active_allunits ;
inclusion_list.amplitude_scale_bar_val = amplitude_scale_bar_val;
inclusion_list.peak_to_peak_amp_micro_volts = peak_to_peak_amp_micro_volts;
% inclusion_list.NlxHeader_Threshold = NlxHeader_Threshold;
inclusion_list.experiment_parameters = experiment_parameters_allunits;
inclusion_list.cell_parameters = cell_parameters_allunits;
inclusion_list.nominal_depth_all_units = nominal_depth_all_units;
inclusion_list.session_name_list = session_name_list_all;
inclusion_list.behavior_sessions = behavior_sessions_list;
inclusion_list.session_start_end_time=session_start_end_time;
inclusion_list.Sync_by_LED =Sync_by_LED;
% inclusion_list.initial_smoothing_factor= P(1,1).LEDs.smoot;

save(filename_out, 'inclusion_list');


