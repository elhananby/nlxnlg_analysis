function NLG_create_incl_list_cell_struct_homing_vector


%---------------------------------------------
% Arseny Finkelstein
%---------------------------------------------

% ========= Stores the individual data for every cell in the inclusion list into a large structure of variables, later (end of the script) to be saved into a mat-file: ========

% ------------- General Parameters: ----------------

p.path_dataout= 'D:\Ayelet\Data\Data_Nlg_Proc\Homing_vector\';

% --------------------------------------------------


% ======= Extract the list of Active Neurons: =======
load([p.path_dataout '\inclusion_list_of_neurons.mat'])  ; % Loading Inclusion-List of Neurons
% This mat-file is where I save all the important data, including:
% * inclusion_list.file_list_allunits = Ntt file names of all my spike-sorted units.
% * inclusion_list.file_list_associated_VT_files_allunits = associated VT files (mat-files).
% * inclusion_list.number_of_spikes_in_every_session_allunits = Number of spikes in every session
% * inclusion_list.is_unit_is_pyramidal_is_active_allunits = Assignments of neurons:
%       [0 1 0] = Not an acceptable unit (e.g. too few spikes overall, or funky waveform).
%       [1 1 0] = An acceptable cell, but is NOT active during behavioral sessions.
%       [1 1 1] = An acceptable AND active cell. ONLY FOR THESE CELLS WILL I COMPUTE PLACE FIELDS AND SPATIAL-VIEW FIELDS.
% * inclusion_list.x_bin_size_pixels_all_units = x-Bin size for computing place fields
% * inclusion_list.y_bin_size_pixels_all_units = y-Bin size for computing place fields
% * inclusion_list.inclusion_list_criteria_comments = the criteria I used for assigning neurons as [1 1 0] etc.

dir_save_data=[p.path_dataout '\Homing_vector_mat\bat_' num2str(inclusion_list.bat(1)),'\']; %Here I will save the data
dir_save_figs = [p.path_dataout '\Homing_vector_flight\Homing_vector_plots_half_sessions_flight/bat_' num2str(inclusion_list.bat(1)),'/']; % Here I will save the figures
Nlg_TTL_delay_ms=35; %the TTL are registered 35 ms after the spikes




idx_acceptable_and_active_cells = ... % List of Active Neurons
    find( inclusion_list.is_unit_is_pyramidal_is_active_allunits(:,3) == 1 );

% ======== Loop over Acceptable+Active cells: ========

for ii_cell = 1:length( idx_acceptable_and_active_cells ), % Loop over Acceptable+Active cells
    
    idx_inclusion_list = idx_acceptable_and_active_cells( ii_cell );
    filename_spike = inclusion_list.file_list_allunits{ idx_inclusion_list };
    filename_VT = inclusion_list.file_list_associated_VT_files_allunits{ idx_inclusion_list };
    session_names = inclusion_list.session_name_list {idx_inclusion_list};
    behavior_sessions = inclusion_list.behavior_sessions{idx_inclusion_list};
    
    %Information from inclusion list
    cell_data.CELL_NUMBER_flight=inclusion_list.cell_parameters{idx_inclusion_list}.cell_num;
    cell_data.BAT_NUMBER=inclusion_list.bat(idx_inclusion_list);
    cell_data.BAT_NAME=inclusion_list.bat_name{idx_inclusion_list};
    cell_data.DAY_flight=inclusion_list.day(idx_inclusion_list);
    cell_data.NOMINAL_DEPTH_MICRON=inclusion_list.nominal_depth_all_units(idx_inclusion_list);
    
    cell_data.PARAMETERS_recording_flight.TT=inclusion_list.TT(idx_inclusion_list);
    cell_data.PARAMETERS_recording_flight.cell_id_on_the_tetrode=inclusion_list.cell_id(idx_inclusion_list);
    cell_data.PARAMETERS_recording_flight.reference_channel=inclusion_list.experiment_parameters{idx_inclusion_list}.reference_channel;
    cell_data.PARAMETERS_recording_flight.is_unit_is_pyramidal_is_active=inclusion_list.is_unit_is_pyramidal_is_active_allunits(idx_inclusion_list,:);
    cell_data.PARAMETERS_recording_flight.peak_to_peak_amp_micro_volts=inclusion_list.peak_to_peak_amp_micro_volts(idx_inclusion_list);
    cell_data.PARAMETERS_recording_flight.r_threshold_for_artifact_cleaning=inclusion_list.experiment_parameters{idx_inclusion_list}.r_threshold;
    
    
    cell_data.ALL_SESSIONS_INFO_flight.pre_analysis_comments=inclusion_list.experiment_parameters{idx_inclusion_list}.comments;
    cell_data.ALL_SESSIONS_INFO_flight.post_analysis_comments=inclusion_list.cell_parameters{idx_inclusion_list}.comments;
    cell_data.ALL_SESSIONS_INFO_flight.throw_away_times=inclusion_list.experiment_parameters{idx_inclusion_list}.throw_away_times;
    cell_data.ALL_SESSIONS_INFO_flight.nsessions=inclusion_list.experiment_parameters{idx_inclusion_list}.nsessions;
    cell_data.ALL_SESSIONS_INFO_flight.sessions_used_for_sorting=inclusion_list.experiment_parameters{idx_inclusion_list}.use_for_sorting;
    cell_data.ALL_SESSIONS_INFO_flight.sessions_used_for_analysis=inclusion_list.experiment_parameters{idx_inclusion_list}.use_for_analysis;
    cell_data.ALL_SESSIONS_INFO_flight.number_of_spikes_in_every_session = inclusion_list.number_of_spikes_in_every_session_allunits( idx_inclusion_list, : ) ;
    
    cell_data.DATA_BY_SESSION = [];
    %% Ayelet: enter for each cell the goals position in the data struct:
    
    excel_sheet=  'inclusion_list_flight.xls';
    sheet=4; %the goal position is the the 4th sheet in the excel
 
    [num,txt,raw] = xlsread(excel_sheet,sheet);
    numeric_fields={'bat','day','central_goal_position','hidden_goal_position_from_event'};
field_names = txt(1,:);
rows = 1:9999;
sheet = raw(2:end,:);
nrows = size(sheet,1);
rows(rows > nrows) = [];
sheet = sheet(rows,:);

% perform eval on numeric fields (surrounded by '[' and ']')

nfields = length(field_names);
field_names_numeric = regexprep(field_names,'\d','#');
numeric_inds = find(ismember(field_names_numeric,numeric_fields));
character_inds = setdiff(1:nfields,numeric_inds);
for i = numeric_inds
    for j = 1:size(sheet,1)
        if ~isnumeric(sheet{j,i})
            sheet{j,i} = eval(['[ ' sheet{j,i} ' ]']);
        end
    end
end
for i = character_inds % turn NANs in charater fields to empty strings
    for j = 1:size(sheet,1)
        if isnan(sheet{j,i})
            sheet{j,i} = '';
        end
    end
end

% convert to struct


   goal_position_struct = cell2struct(raw,field_names,2);
    %find the goal position for the relevand day:
    
    day_line_num=find(strcmp(raw,num2str(cell_data.DAY_flight)))-(size(raw,1));
    
    
    cell_data.ALL_SESSIONS_INFO_flight.central_goal_position=str2num(goal_position_struct(day_line_num,1).central_goal_position);
    cell_data.ALL_SESSIONS_INFO_flight.hidden_goal_position_from_event=str2num(goal_position_struct(day_line_num,1).hidden_goal_position_from_event);
    


   
%%
    
    for ii_session = 1:length(session_names) % loop over behavioral sessions for analysis
        session_name = session_names{ii_session};
        session_number=ii_session;
        
        load(filename_VT{ii_session}); % Load VT (video) data for session
        
        % ==== Read the Ntt files, and extract the data that occureed within the session: ====
        FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
        ExtractionMode = 1 ; % Mode 1 = "Extract All"
        ExtractionModeArray = [] ; % Will read all the data
        
        [Timestamps, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
            Nlx2MatSpike( filename_spike, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
        
        timestamps_session_limits=inclusion_list.session_start_end_time{ii_cell, ii_session};
        
        idx_within_session = find( Timestamps >= timestamps_session_limits(1) & ...
            Timestamps <= timestamps_session_limits(2) ); % Find spikes that occured within this session
        
        SPIKES_Timestamps = Timestamps( idx_within_session );
        
%         %% Substracting 35 msec due to Neurologger delay in registering the TTL pulses
%         if inclusion_list.Sync_by_LED(ii_cell)==0 %i.e it was synchronized by TTLs
%             SPIKES_Timestamps=SPIKES_Timestamps-(Nlg_TTL_delay_ms*1000);
%         end;
        
%         COMPUTE HEADING DIRECTION
%         [VT_Timestamps_cleaned_FE, SPIKES_Timestamps_cleaned_FE, DATA_BY_SESSION_flight] ...
%             = heading_direction_flight (VT, VT_Parameters, SPIKES_Timestamps, timestamps_session_limits, inclusion_list, idx_inclusion_list,dir_save_figs);
%         
        % Storing the Data per cell/session

%         cell_data.DATA_BY_SESSION{session_number}=DATA_BY_SESSION_flight;
        
        cell_data.DATA_BY_SESSION{session_number}.session_name = session_names(ii_session) ;
        %cell_data.DATA_BY_SESSION{session_number}.total_number_of_spikes = length(SPIKES_Timestamps_cleaned_FE) ;
        
        %cell_data.DATA_BY_SESSION{session_number}.position.Timestamps = VT_Timestamps_cleaned_FE;
        cell_data.DATA_BY_SESSION{session_number}.spikes.Timestamps = SPIKES_Timestamps;
         
%         % Although we don't measure this parameters, we just  fill in this to be compatible with the crawling code        
%         cell_data.DATA_BY_SESSION{session_number}.number_of_spikes_with_good_tracking=cell_data.DATA_BY_SESSION{session_number}.total_number_of_spikes;
%         cell_data.DATA_BY_SESSION{session_number}.position.roll_angle=cell_data.DATA_BY_SESSION{session_number}.position.pitch_angle;
%         cell_data.DATA_BY_SESSION{session_number}.spikes.roll_angle_at_spike=0*cell_data.DATA_BY_SESSION{session_number}.spikes.pitch_angle_at_spike;
%         cell_data.DATA_BY_SESSION{session_number}.spikes.roll_angle_at_spike_by_closest_frame=cell_data.DATA_BY_SESSION{session_number}.spikes.pitch_angle_at_spike_by_closest_frame;
%         
        % We also save this parameters (for flight session only) needed for
        % 3D place cells calculations
        cell_data.DATA_BY_SESSION{session_number}.VT=VT;
        cell_data.DATA_BY_SESSION{session_number}.VT_Parameters=VT_Parameters;
        cell_data.DATA_BY_SESSION{session_number}.timestamps_session_limits=timestamps_session_limits;
        
        cell_data.COMMENT_flight = ...
            'Data saved by the script:  D:\Arseny\matlab_code\Projects\3D_HD_flight\NLG_create_incl_list_cell_struct.m AND heading_direction_flight.m' ;
%% Ayelet- add the 2nd session timestamps:
if ii_session==2
% find the string of start/end session in the events strings:         
find_end_session_1_event=strfind(inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventStrings,'End session 1');
find_start_session_2_event=strfind(inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventStrings,'Start session 2');

end_session_1_event=find(~cellfun(@isempty,find_end_session_1_event));
start_session_2_event=find(~cellfun(@isempty,find_start_session_2_event));
%find the event timestamp:
end_session_1_timestamp=inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventTimestamps(end_session_1_event);
start_session_2_timestamp=inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventTimestamps(start_session_2_event);
   
%Save the session timestamps into the cell_data struct:
cell_data.DATA_BY_SESSION{session_number}.end_session_1_timestamp=end_session_1_timestamp(end);        
cell_data.DATA_BY_SESSION{session_number}.start_session_2_timestamp=start_session_2_timestamp(end);    
 find_start_session_3_event=strfind(inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventStrings,'Start session 3');

% find for session 3:
if ~isempty(find(~cellfun(@isempty,find_start_session_3_event)))
find_end_session_2_event=strfind(inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventStrings,'End session 2');
find_start_session_3_event=strfind(inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventStrings,'Start session 3');

end_session_2_event=find(~cellfun(@isempty,find_end_session_2_event));
start_session_3_event=find(~cellfun(@isempty,find_start_session_3_event));
%find the event timestamp:
end_session_2_timestamp=inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventTimestamps(end_session_2_event);
start_session_3_timestamp=inclusion_list.experiment_parameters{idx_inclusion_list}.Nlg_EventTimestamps(start_session_3_event);
   
%Save the session timestamps into the cell_data struct:
cell_data.DATA_BY_SESSION{session_number}.end_session_2_timestamp=end_session_2_timestamp(end);        
cell_data.DATA_BY_SESSION{session_number}.start_session_3_timestamp=start_session_3_timestamp(end);    
  
end  
end
    end % enumerating over behavior sessions
    
    
    
%     %% ======= Combine the Flight session with the Crawling sessions (if exist) ==========
%     
%     dir_data_HD_in_crawling='D:\Arseny\Data_PrS_ring\Head_Direction\';
%     b1_similar_cell= inclusion_list.cell_parameters{idx_inclusion_list}.b1_similar_cell;
%     b3_similar_cell= inclusion_list.cell_parameters{idx_inclusion_list}.b3_similar_cell;
%     
%     % update the data (if one of the crawling sessions exist)
%     if b1_similar_cell ~=0
%         current_cell=dir(fullfile([dir_data_HD_in_crawling  'HD_cell#' num2str(b1_similar_cell)  '*.mat']));
%         b1 = load([dir_data_HD_in_crawling current_cell.name]);
%         cell_data.DATA_BY_SESSION{1}=b1.DATA_BY_SESSION{1};
%     end
%     if b3_similar_cell ~=0
%         current_cell=dir(fullfile([dir_data_HD_in_crawling  'HD_cell#' num2str(b3_similar_cell)  '*.mat']));
%         b3 = load([dir_data_HD_in_crawling current_cell.name]);
%         cell_data.DATA_BY_SESSION{3}=b3.DATA_BY_SESSION{end};
%     end
%     
%     % update the parameters (if one of the crawling sessions exist)
%     if b1_similar_cell ~=0
%         cell_data.CELL_NUMBER = b1.CELL_NUMBER;
%         cell_data.PARAMETERS_fixed_for_all_units = b1.PARAMETERS_fixed_for_all_units;
%         cell_data.DAY = b1.DAY;
%         cell_data.NEURALYNX_FILE_LOCATION = b1.NEURALYNX_FILE_LOCATION;
%         cell_data.PARAMETERS_recording = b1.PARAMETERS_recording;
%         cell_data.ALL_SESSIONS_INFO = b1.ALL_SESSIONS_INFO;
%         cell_data.PARAMETERS_video_tracking = b1.PARAMETERS_video_tracking;
%         cell_data.COMMENT = b1.COMMENT;
%         
%     elseif b3_similar_cell ~=0
%         cell_data.CELL_NUMBER = b3.CELL_NUMBER;
%         cell_data.PARAMETERS_fixed_for_all_units = b3.PARAMETERS_fixed_for_all_units;
%         cell_data.DAY = b3.DAY;
%         cell_data.NEURALYNX_FILE_LOCATION = b3.NEURALYNX_FILE_LOCATION;
%         cell_data.PARAMETERS_recording = b3.PARAMETERS_recording;
%         cell_data.ALL_SESSIONS_INFO = b3.ALL_SESSIONS_INFO;
%         cell_data.PARAMETERS_video_tracking = b3.PARAMETERS_video_tracking;
%         cell_data.COMMENT = b3.COMMENT;
%     end;
    
    %% ======= Save mat-file with the data: ==========
    
    
    
    if isempty(dir(dir_save_data))
        mkdir(dir_save_data);
    end;
    
    filename= [dir_save_data 'HV_cell#' num2str(inclusion_list.cell_parameters{idx_inclusion_list}.cell_num)...
        '_bat' num2str(inclusion_list.cell_parameters{idx_inclusion_list}.bat)...
        '_TT' num2str(inclusion_list.cell_parameters{idx_inclusion_list}.TT)...
        '_unit' num2str(inclusion_list.cell_parameters{idx_inclusion_list}.cell_id)...
        '_' num2str(inclusion_list.cell_parameters{idx_inclusion_list}.day)...
        '.mat'] ; % I will save the data into this mat-file
    
    save (filename, '-struct', 'cell_data'); %saves only the subfields of cell_data
    
end % enumerating over the cells



% PrS_analyze_angles_flight_merged_sessions([]);










