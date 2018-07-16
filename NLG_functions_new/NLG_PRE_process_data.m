%---------------------------
% Arseny Finkelstein 19/03/2014

% Some parts of the code don't work in versions older then 2012b

% batch program for pre-processing Nlg and Nlx information
% read excel file containing info
%
% input:
% excel_sheet - containing rows of data about experiemental sessions
%               (see for example inclusion_list_sessions.xls)
%               Note that excel sheet is assumed to be in 'text' format.
% rows - which rows to read in Experiments sheet. Default is 1-9999
% cell_rows - which rows to read in Cells sheet. Default is 1-9999
%--------------------------
close all;
clear all;

%% open matlab workers (for parfor)
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 7
else
    disp('matlab workers already open')
end

%%%%% general parameters that do not change from file to file
p_in.Nlx_timestamps=[]; % will be filled automatically later on
p_in.Spike_threshold_uV_units = 50; % Used as threshold for spike detection in Nlg_Flight_detect spikes CSC
p_in.include_negative_threshold=1; % 1= take 'inverted' spikes - with valley more negative than ( - Spike_threshold_uV_units)
p_in.do_comparison_acceptable_shapes=1; % 1= run comparison to the library of acceptable shapes; 0 - skipt this step
p_in.do_coincidence_detection=1; % 1= run coincidence detection; 0 - skip this step
p_in.r_threshold = 0.8;   %  Analog system value is 0.92
p_in.remove_bad_frames=0;  %1= remove bad frames and not detect spikes there.
p_in.path_datain = 'D:\Ayelet\Data\'; % prefix for all data folders
p_in.path_year_bat=[]; % will be filled automatically later on
p_in.path_dataout = 'D:\Ayelet\Data\Data_Nlg_Proc\'; % prefix for output folders
p_in.numeric_fields = {'bat','day','reference_channel','TT','depth', ...
    'active_channels','nsessions','use_for_sorting', 'use_for_analysis', ...
    'use_tetrodes','throw_away_times', ...
    'events_#','time_offsets_#_in_seconds','exact_time_in_microseconds_#' 'Nlx_events', 'Nlx_time_offsets_in_seconds','Nlx_exact_time_in_microseconds','Sync_by_LED',...
    'cell_num', 'cell_id','b1_similar_cell','b3_similar_cell', 'single_unit','pyramidal', 'behaviorally_active'}; % numeric fields in excel
p_in.Nlg_EventStrings=[]; % will be filled automatically later on
p_in.Nlg_EventTimestamps=[]; % will be filled automatically later on
p_in.Nlx_EventStrings=[]; % will be filled automatically later on
p_in.Nlx_EventTimestamps=[]; % will be filled automatically later on
p_in.TTL_timestamps_Nlg=[]; % will be filled automatically later on
p_in.TTL_timestamps_Nlx=[]; % will be filled automatically later on
p_in.polyfit_Nlg2Nlx_microsec=[]; % will be filled automatically later on
p_in.polyfit_Nlx2Nlg_microsec=[]; % will be filled automatically later on
% ------------- General Parameters for video tracker: ----------------


%% Read excel sheet and create parameter structure from it

if ~exist('excel_sheet') || isempty(excel_sheet) % defaults
    excel_sheet = 'inclusion_list_flight.xls';
end
if ~exist('rows') || isempty(rows)
    rows = 1:9999;
end
if ~exist('cell_rows') || isempty(cell_rows)
    cell_rows = 1:9999;
end

P = PRE_read_excel_sheet(excel_sheet,'Experiments',rows,p_in.numeric_fields,p_in);


if length(P) == 0 %if no records are entered
    return
end


%% Extract events for each day

for ii_rec = 1:length(P) % number of records (days) in excel
    
    p = P(ii_rec);
    
    %Generating path for the data folder
    day_string= num2str(p.day);
    year=num2str(P(1).day);%fix to bug change in years!
    %DEBUG- change it after artifact test:
    p.path_year_bat=['yr' year(1:4) '_bat' num2str(p.bat) '_' P(1, 1).bat_name '_Nlg\'];
    %p.path_year_bat=['yr' day_string(1:4) '_bat' num2str(p.bat) '_' p.bat_name '_Nlg_test_artifacts\'];

    disp('========================================================');
    disp(' ');
    disp(['Extracting day' num2str(ii_rec) ' ' day_string]);
    disp(' ');
    
    
    %----------------------------------------------------------------------
    % 1) Extracting ans Show Events and Timestamps for both Nlg and Nlx data
    %----------------------------------------------------------------------
    
    
    %Neurologger: Extracting Events, Event Timestamps, and TTL times from the event log of the Nlg
    %     file_name_excel_Event_Nlg=  [ p.path_datain p.path_year_bat p.path_day_dir '\RAW_DATA\EVENTLOG.xlsx'];
    file_name_Event_Nlg=  [ p.path_datain p.path_year_bat p.path_day_dir '\nlx\'];
    [p.Nlg_EventStrings, p.Nlg_EventTimestamps, p.TTL_timestamps_Nlg] = NLG_PRE_read_excel_events_and_TTL(file_name_Event_Nlg);
    
    
    disp('Neurologger EVENT LIST:');
    disp('===========');
    p.Nlg_EventStrings % displaying event list
    
    %Neuralynx: : Extracting Events
    disp('Neuralynx EVENT LIST:');
    disp('===========');
    [p.Nlx_EventStrings p.Nlx_EventTimestamps] =NLG_PRE_show_Nlx_events(p); % display event list
    
    
    %----------------------------------------------------------------------
    % 2) Check that all info is filled in in excel
    %----------------------------------------------------------------------
    
    %1) check that all Nlg sessions are defined in excel
    if ~isfield(p,'S') || isempty(p.S)
        disp('please complete Nlg session names info in excel and run again');
        return;
    end
    %2) check that all Nlg sessions events are defined in excel
    nsessions = length(p.S);
    for nses = 1:nsessions
        s = p.S(nses);
        if isempty(s.events) || all(isnan(s.events))
            disp('please complete Nlg session events info in excel and run again');
            return;
        end;
    end;
    %3) check that Nlx information is defined in excel
    if isempty(p.Nlx_session) || isempty(p.Nlx_session)
        disp('please complete Nlx info excel and run again');
        return;
    end;
    
    %----------------------------------------------------------------------
    % 3) For each session defined in excel get the following parameters:
    %----------------------------------------------------------------------
    
    % start_behavior  - pointer to start of behavior in event list
    % end_behavior    - pointer to end of behavior in event list
    % event_list      - event list itself
    % time_stamps     - time stamps of each event in event list
    % time_stamps offsets - if entered manually by the user
    
    % Important:
    %Nlx - timestamps in microsec
    %Nlg - timestamps also in microseconds (!)
    
    p = NLG_PRE_get_session_times_Nlg_and_Nlx(p);
    
    
    %----------------------------------------------------------------------
    % 9) Update the structure (inclusion list)
    %----------------------------------------------------------------------
    P(ii_rec) = p;
    
end % loops on excel records for display of events  only


%% Extract Timestamps, TTLs, Video and Spikes, for each day
% 
for ii_rec = 1:length(P) % number of records (days) in excel
    p = P(ii_rec);
    
    %Generating path for the data folder
    day_string= num2str(p.day);    
    year=num2str(P(1).day); %fix to bug change in years!
    p.path_year_bat=['yr' year(1:4) '_bat' num2str(p.bat) '_' P(1, 1).bat_name '_Nlg\'];
     % p.path_year_bat=['yr' day_string(1:4) '_bat' num2str(p.bat) '_' p.bat_name '_Nlg_test_artifacts\'];

    
    disp('========================================================');
    disp(' ');
    disp(['Extracting day' num2str(ii_rec) ' ' day_string]);
    disp(' ');
    
    %----------------------------------------------------------------------
    % 4) Extract Nlg CSC data, filter it, and divide into 1-min chunks
    %----------------------------------------------------------------------
    Nlg_Flight_extract_and_filter_CSC_parfor(p); %runs on all channels including the 'bad' ones
    
    %----------------------------------------------------------------------
    % 5) Clean CSC Artifacts epochs
    %----------------------------------------------------------------------
    Nlg_clean_artifacts_CSC_parfor(p); %runs on all channels including the 'bad' ones
    
    %----------------------------------------------------------------------
    % 6) Sync Nlx2Nlg by either TTLs events or LEDs on/off times
    %----------------------------------------------------------------------
    folder_name=[p.path_datain p.path_year_bat p.path_day_dir];
    
    if p.Sync_by_LED==0 % We are synchronizing by TTLs
        %FindTTL timestamps detected in the Nlx clock (in microsec)
        % Extracting TTL onset times from Nlx, by crosscorrelation-threshold
        % detection of the continiously sampled TTL channel ('Audio2')
        file_name_Nlx_TTL=[p.path_datain p.path_year_bat p.path_day_dir '\Nlx_VT_and_TTL\Audio2.ncs'];
             %   [error_sync]= NLG_PRE_sync_Nlx2Nlg_by_TTL_2_0(file_name_Nlx_TTL,p,folder_name)

        [p.TTL_timestamps_Nlx  p.polyfit_Nlg2Nlx_microsec p.polyfit_Nlx2Nlg_microsec error_sync]= NLG_PRE_sync_Nlx2Nlg_by_TTL (file_name_Nlx_TTL,folder_name,p);
        close all;
        
        %If there was a problem with synchronization we will skip processing
        %this day
        if error_sync==1
            continue
        end;
        
    elseif p.Sync_by_LED==1 % We are synchronizing by LEDs on/off times
        %FOR DBG - to find the delay between the emission of TTL pulse and the registration of its timespent on Nlg
        %file_name_Nlx_TTL=[p.path_datain p.path_year_bat p.path_day_dir '\Nlx_VT_and_TTL\Audio2.ncs'];
        %NLG_DBG_find_TTL_delay_by_LED(file_name_Nlx_TTL,folder_name, p);
        [p.polyfit_Nlg2Nlx_microsec p.polyfit_Nlx2Nlg_microsec ]= NLG_PRE_sync_Nlx2Nlg_by_LED (folder_name, p);
        close all;
    end
    
    %----------------------------------------------------------------------
    % 7)Extract the Video and compute the trajectory in 3D
    % (by  synchronizing the Nlx Video-tracker into Nlg timestampts
    %----------------------------------------------------------------------
    
    Nlg_extract_video_flight(p); %CSC artifacts are removed from the coverage data based on 'good' channels only
    close all;
    %----------------------------------------------------------------------
    % 7)Detect spikes from CSC
    %----------------------------------------------------------------------
    
   Nlg_Flight_detect_spikes_CSC_parfor(p); %CSC artifacts are removed from the neural data based on 'good' channels only
    
    %----------------------------------------------------------------------
    % 9) Update the structure (inclusion list)
    %----------------------------------------------------------------------
    P(ii_rec) = p;
    
end % loops on excel records for data (video and spikes extraction)


% check whether user has filled cells in excel sheet 

% create cell inclusion list

C = PRE_read_excel_sheet(excel_sheet,'Cells',cell_rows,p_in.numeric_fields,[]);

if isempty(C)
    return;
end

NLG_PRE_create_inclusion_list_of_neurons(P,C);
NLG_create_incl_list_cell_struct_homing_vector;