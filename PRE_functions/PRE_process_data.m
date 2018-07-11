function PRE_process_data(varargin)
%
% batch program for pre-processing neuralynx information
% read excel file containing info
%
% input:
% excel_sheet - containing rows of data about experiemental sessions
%               (see for example inclusion_list_sessions.xls)
%               Note that excel sheet is assumed to be in 'text' format.
% rows - which rows to read in Experiments sheet. Default is 1-9999#

% Example:
% PRE_process_data('D:\inclusion_list_sessions.xls',[],1) processes only cells in rows 1 and 2 in the
%                           excel sheet
%
dbstop if error;

%% input processing
switch nargin
    case 0 % if no variables
        excel_sheet = 'D:\Scripts\nlxnlg_analysis\inclusion_list_sessions.xlsx';
        T_temp = readtable(excel_sheet, 'Sheet', 1);
        rows = 1:height(T_temp);
        
    case 1 % if one variable (only excel sheet)
        if ~isempty(varargin{1})
            excel_sheet = varargin{1};
        else
            excel_sheet = 'D:\Scripts\nlxnlg_analysis\inclusion_list_sessions.xlsx';
        end
        T_temp = readtable(excel_sheet, 'Sheet', 1);
        rows = 1:height(T_temp);
        
    case 2 % if two variables (
        if ~isempty(varargin{1})
            excel_sheet = varargin{1};
        else
            excel_sheet = 'D:\Scripts\nlxnlg_analysis\inclusion_list_sessions.xlsx';
        end
        
        if ~isempty(varargin{2})
            rows = varargin{2};
        else
            T_temp = readtable(excel_sheet, 'Sheet', 1);
            rows = 1:height(T_temp);
        end
        
end

%% setting parameters
p_in.r_threshold = 0.9;   % the higher this value, less spikes are accepted
p_in.path_datain = '\\GutfreundNLX\CheetahData'; % prefix for all data folders
p_in.path_dataout ='D:\experiment_data'; % prefix for output folders
p_in.cell_dir = 'cells';  % directory of cells database (under Data_p)
p_in.cell_prefix = 'quail'; % prefix for cell names in database

% ------------- General Parameters for video tracker: ----------------
p_in.VT_Resolution = [720 576]; % Number of pixels in VT image: [X Y]
p_in.num_lightsticks_Arena_Marking = 2 ; % Number of light markers (LED's in arena corners) I used for Arena-Marking ( = No.of Targets to read)
p_in.boxSize = 96.5 ; % The distance between the East to West walls of the arena (in cm)
p_in.smoothing_parameter_for_velocity = 10^-5 ; % STRONG smoothing for Velocity computation (csaps.m/fnder.m)

% fields considered as numeric when reading excel

numeric_fields = {'animal', 'day','experiment','reference_channel','TT','depth', ...
    'active_channels','nsessions','use_for_sorting', 'use_for_analysis', ...
    'use_tetrodes', 'NW_calib', 'NE_calib', 'SE_calib', 'SW_calib', ...
    'events_#','time_offsets_#', 'cells', 'cell_number', ...
    'cell_id', 'single_unit','pyramidal', 'behaviorally_active',...
    'throw_away_times','arena_width_east_to_west'}; % numeric fields in excel %add 'arena_#' to the list

P = PRE_read_excel_sheet(excel_sheet, 'Experiments', rows, numeric_fields,p_in);
nrecs = length(P); % number of records in excel   %replace length(P) with 1 if not more than one row

if nrecs == 0
    return
end

% check that all events are defined in excel and loop on excel records
for nrec = 1:nrecs
    
    p = P(nrec);
    
    % disaply calibration event list
    if isempty(p.NW_calib)
        PRE_show_events(p, 'calibration');
    end
    
    % display event list
    if ~exist('p.S(1).events', 'var')
        PRE_show_events(p, 'recording');
    end
    
    if ~isfield(p,'S') || isempty(p.S)
        disp('please complete session info in excel and run again');
        return;
    end
    
    nsessions = length(p.S);
    
    for nses = 1:nsessions
        s = p.S(nses);
        if isempty(s.events) || all(isnan(s.events))
            disp('please complete session info in excel and run again');
            return;
        end
    end
    
    p = PRE_get_session_times(p);
    
    if strcmp(p.nlgnlx, 'nlg')
        p.datadir_out = sprintf('nlg_animal%d_Day%d_%d\\', p.animal, p.day, p.experiment);
        datadir_out = fullfile(p.path_dataout, p.datadir_out);
        p = NLG_PRE_process_data(p);
    else
        p.datadir_out = sprintf('animal%d_Day%d_%d\\', p.animal, p.day, p.experiment);
        datadir_out = fullfile(p.path_dataout, p.datadir_out);
    end    

    %% clean spike artifacts
    if ~strcmp(p.nlgnlx, 'nlg') % if its a nlg recording, we already did spike cleaning
        PRE_clean_artifacts(p, datadir_out);
        p.nlg = [];
    end
    
    %% extract data from video
    p = PRE_extract_video_w_reflection_fix(p);

    P_new(nrec) = p;
    
end % loops on excel records

P = P_new;

%% pause
fprintf('Finished spike detection and video analysis,\nPerform spike sorting and press any key to continue\n');
% pause();

%% read cells from sorted ntt file
cells = PRE_create_cells_list(excel_sheet, P);

%% find cells to load
cellsToProcess = subdir('D:\experiment_data\Cells\*.mat');
idx = 1;
cells = {};
printLine = 0;

for iiCell = 1:length(cellsToProcess)
    % display progress
    fprintf(repmat('\b', 1, printLine));
    printLine = fprintf('%.2f', iiCell*100/length(cellsToProcess));
    
    % load cell
    load(cellsToProcess(iiCell).name);
    
    % calculate occupancy percentage
    binEdges = linspace(0, 100, 30);
    occupancyMap = histcounts2(vt.posx_c, vt.posy_c, binEdges, binEdges);
    occupancyPercentage = numel(find(occupancyMap))./numel(occupancyMap);
    
    % calculate speed score
    [~, ~, ~, speedScore] = calculate_speed_map(c, vt);
    
    % keep only cells where the session >= 10 minutes and more than 300
    % spikes and occupancy > .75
    if (p.S(metaData.session).end_time - p.S(metaData.session).start_time)*1e-6/60 >= 10 &&...
            length(c.timestamps) >= 300 && abs(speedScore >= .2)
        cells{idx} = cellsToProcess(iiCell).name;
        idx = idx + 1;
    end
end

%% do a very basic analysis on the data
basic_analysis(cells);

fprintf('finished complete run of PRE_process_data');
end