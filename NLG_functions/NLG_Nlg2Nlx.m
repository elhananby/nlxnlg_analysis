function NLG_Nlg2Nlx(p)

header_file = 'D:\Scripts\nlx_analysis\NLG_functions\header.txt';
Nlg_InDir = fullfile(p.path_datain, p.data_dir, 'nlg_data');
Nlx_OutDir = fullfile(p.path_dataout, p.datadir_out, 'nlx_data');

if exist(Nlx_OutDir, 'dir') && ~isempty(subdir(fullfile(Nlx_OutDir, '*.dat')))
    fprintf('NLG files already processed - skipping.\n');
    return;
end

num_channels = 16;
data_cnl_ind = [0:15]; % note that this numbering system is of the neurologger which means channels 0-15
DATA_file_prefix = 'NEUR';
zero_DC_level_bit = 2048;
is_invert_data = true;
is_remove_DC = true;
is_remove_flash_write_artifact = true;
use_clock_diff_correction = true;
use_post_rec_ref_channel = false; % use this if you recorded with GND as ref channel
%post_rec_ref_channel = 10; % If the above is true - choose the channel you want to substract from all the other channels

% chanles number are from 1 to 16
%% read EVENT file
event_file_name_xlsx = fullfile(Nlg_InDir, 'EVENTLOG.CSV');
[NUM, TXT, RAW]=xlsread(event_file_name_xlsx);

% extract recording details from event file header
file_header_lines = TXT(1:3,1);
[splitstr] = regexp(file_header_lines{2}, '[;]+', 'split'); % 2nd header row
firmware_ver = regexp(splitstr{1}, '\d*[.]\d*','match');
serial_number = regexp(splitstr{2}, '\d*','match');
time = regexp(splitstr{3}, '\d*:\d*:\d*','match');
date = regexp(splitstr{4}, '\d*/\d*/\d*','match');
ADC_period_usec = regexp(splitstr{5}, '\d*[.]\d*','match');

if isempty(ADC_period_usec)
    ADC_period_usec = regexp(splitstr{5}, '(\d*)', 'match');
end
[splitstr] = regexp(file_header_lines{3}, '[;]+', 'split'); % 3rd header row
ADC_resolution_uVolt = regexp(splitstr{1}, '\d*[.]\d*','match');

ADC_SAMPLE_PERIOD = str2num(cell2mat(ADC_period_usec))/num_channels*1e-6;
fs =  1/(ADC_SAMPLE_PERIOD * num_channels);
uVolt_per_bit = str2num(cell2mat(ADC_resolution_uVolt));
block_period_time_usec = (512/fs) * 1e6;
% block_period_time_usec_2 = 512 * (ADC_SAMPLE_PERIOD * num_channels) * 1e6;


save(fullfile(Nlx_OutDir, 'params.mat'));
%% extract events details
events_IX = NUM(1:end,1);
events_TS = NUM(1:end,3).*1e3;
events_TS_source = TXT(5:end,4);
events_type = TXT(5:end,5);
events_details = TXT(5:end,6);

%% join '...Continued' events
continued_event_lines_IX = find(isnan(events_IX))';
for ii_continued_event_line = continued_event_lines_IX
    % take the last line with a valid number in the event index column as
    % the event index
    last_valid_line = find( ~isnan(events_IX(1:ii_continued_event_line)), 1, 'last');
    events_details{last_valid_line} = [events_details{last_valid_line} events_details{ii_continued_event_line}];
end

events_IX(continued_event_lines_IX) = [];
events_TS(continued_event_lines_IX) = [];
events_TS_source(continued_event_lines_IX) = [];
events_type(continued_event_lines_IX) = [];
events_details(continued_event_lines_IX) = [];

%% join event type with event details to a single string
EventStrings = {};
for ii_event = 1:length(events_type)
    event_str = [events_type{ii_event} ' - ' events_details{ii_event}];
    EventStrings{ii_event} = event_str;
end

%% Apply clock difference correction (Transceiver vs. logger time)
if use_clock_diff_correction
    
    PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
    CD_str = 'CD=';
    CD_values = [];
    CD_timestamps = [];
    
    for ii_event = 1:length(PC_gen_events_IX)
        curr_event_details = events_details{PC_gen_events_IX(ii_event)};
        CD_str_pos = strfind(curr_event_details, CD_str);
        
        if isempty(CD_str_pos)
            continue;
        end
        
        CD_values_str_interval = (CD_str_pos:CD_str_pos+4) + length(CD_str);
        CD_values(end+1) = str2num( curr_event_details(CD_values_str_interval) );
        CD_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
    end
    
    CD_values = CD_values.*1e6; % change from sec to usec
    CD_timestamps(find(CD_values>(mean(CD_values)+2*std(CD_values))))=[];
    CD_values(find(CD_values>(mean(CD_values)+2*std(CD_values))))=[];
    
    % TODO: we need to check what is the meaning of CD (clock difference). Logger-Tx
    % or Tx-Logger?
    CD_event_ts__logger_time = CD_timestamps-CD_timestamps(1);        % removed 1st timestamp to balance the fit data closer to zero (otheriwse it arise a warning message)
    CD_event_ts__Tx_time = CD_event_ts__logger_time - CD_values;
    transceiver_2_logger_time_fit = polyfit(CD_event_ts__Tx_time, CD_event_ts__logger_time , 1);
    ts_source_transceiver_events_IX = find(strcmp('Transceiver', events_TS_source));
    events_TS(ts_source_transceiver_events_IX) = polyval(transceiver_2_logger_time_fit , events_TS(ts_source_transceiver_events_IX) - CD_timestamps(1)) + CD_timestamps(1);
    for ii_event = 1:length(ts_source_transceiver_events_IX)
        events_TS_source{ts_source_transceiver_events_IX(ii_event)} = 'Logger';
    end
    
end

%% write event file in Nlx format
if ~exist(Nlx_OutDir, 'dir')
    mkdir(Nlx_OutDir);
end

NlxEventFile = fullfile(Nlx_OutDir, 'EVENTS.nev');
Mat2NlxEV(NlxEventFile, 0, 1, [], [1 0 0 0 1 0] , events_TS', EventStrings');

%% create sperate event file for each event category
event_type_list = unique(events_type);
for ii_event_type = 1:length( event_type_list )
    event_type_string = event_type_list{ii_event_type};
    event_type_events_IX = find(strcmp( event_type_string , events_type));
    Nlx_event_type_file = fullfile(Nlx_OutDir, ['EVENTS__' event_type_string '.nev']);
    Mat2NlxEV(Nlx_event_type_file , 0, 1, [], [1 0 0 0 1 0] , events_TS(event_type_events_IX)', EventStrings(event_type_events_IX)');
end

%% create battery discharge plot
PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
mode_changed_events_IX = find(strcmp('Mode change', events_type));
BV_str = 'BV=';
BV_values = [];
BV_timestamps = [];
for ii_event = 1:length(PC_gen_events_IX)
    curr_event_details = events_details{PC_gen_events_IX(ii_event)};
    BV_str_pos = strfind(curr_event_details, BV_str);
    if isempty(BV_str_pos)
        continue;
    end
    
    BV_values_str_interval = (BV_str_pos:BV_str_pos+4) + length(BV_str);
    BV_values(end+1) = str2num( curr_event_details(BV_values_str_interval) );
    BV_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
end

%% create Nlx data files - continuous sampling files called (.ncs)
% identify  'File started' events and take timestamps
FileStarted_IX = strcmp('File started', events_type);
FileStarted_TS = events_TS(FileStarted_IX);
FileStarted_details = events_details(FileStarted_IX);

% read header template
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
% change fields that are recording specific
header = strrep(header,'-SamplingFrequency', ['-SamplingFrequency ' char(vpa(fs))]);

filesToOpen = subdir(fullfile(Nlg_InDir, '*.DAT'));
Timestamps = [];
fprintf('Writing ');
for ii_file_start_entry = 1:length(FileStarted_TS)
    
    temp = FileStarted_details{ii_file_start_entry};
    file_name = fullfile(filesToOpen(ii_file_start_entry).name);
    channelFid = fopen(file_name);
    filedata = fread(channelFid, '*uint16');
    fclose(channelFid);
    filedata = double(filedata);
    data = reshape(filedata, num_channels, length(filedata)/num_channels);
    
    % remove DC
    if is_remove_DC
        for ii_channel = 1:size(data,1)
            %         data(ii_channel,:) = data(ii_channel,:) - mean(data(ii_channel,:));  % TAMIR: we prefer not to use this since data could be biased around real zero voltage (artifacts for example...)
            data(ii_channel,:) = data(ii_channel,:) - zero_DC_level_bit;         % There should be a theoretical constant representing the real zero voltage, but because the INTAN have a DC dependency on freq it is not perfect
        end
    end
    
    count=0;
    
    % remove repetative flash write artifact
    if is_remove_flash_write_artifact
        weak_artifacts_block_IX = [1:256 513:768 1025:1280 1537:1793];
        strong_artifacts_block_IX = setdiff(1:2048, weak_artifacts_block_IX);
        trim_prc = 5;
        for ii_channel = 1:size(data,1)
            data_channel_blocks = reshape( data(ii_channel,:), 256, [] );
            flash_write_artifact_shape_weak = trimmean(data_channel_blocks(:,weak_artifacts_block_IX), trim_prc , 'round', 2);
            flash_write_artifact_shape_strong = trimmean(data_channel_blocks(:,strong_artifacts_block_IX), trim_prc , 'round', 2);
            data_channel_blocks_corrected = zeros(size(data_channel_blocks));
            
            prc=10:10:90;
            flash_write_artifact_shape_weak_prc=[];
            flash_write_artifact_shape_weak_prc=[];
            for prc_i=1:length(prc)
                
                flash_write_artifact_shape_weak_prc(prc_i,:)=prctile(data_channel_blocks(:,weak_artifacts_block_IX)',prc(prc_i))';
                %flash_write_artifact_shape_weak_prc(prc_i,:)=flash_write_artifact_shape_weak_prc(prc_i,:)-mean(flash_write_artifact_shape_weak(200:end));
                flash_write_artifact_shape_strong_prc(prc_i,:)=prctile(data_channel_blocks(:,strong_artifacts_block_IX)',prc(prc_i))';
                %flash_write_artifact_shape_strong_prc(prc_i,:)=flash_write_artifact_shape_strong_prc(prc_i,:)-mean(flash_write_artifact_shape_strong(200:end));
            end
                   
            for ii_block = weak_artifacts_block_IX
                data_channel_blocks_corrected(:,ii_block) = data_channel_blocks(:,ii_block) - flash_write_artifact_shape_weak;
                
            end
            
            for ii_block = strong_artifacts_block_IX
                data_channel_blocks_corrected(:,ii_block) = data_channel_blocks(:,ii_block) - flash_write_artifact_shape_strong;
            end
            data(ii_channel,:) = reshape(data_channel_blocks_corrected,1,[]);  
        end 
    end
    
    % if we recorded with GND as ref. channel we want to have a
    % "post-recording ref. channel"
    if use_post_rec_ref_channel
        gain_relative_2_ref_channel = 1.*ones(1,16);
        for ii_channel = 1:size(data,1)
            if ii_channel == post_rec_ref_channel
                continue;
            end
            data(ii_channel,:) = data(ii_channel,:) - gain_relative_2_ref_channel(ii_channel).*data(post_rec_ref_channel,:);
   
        end
    end
    
    % change to uVolt units
    data = data.*uVolt_per_bit;

    if is_invert_data
        % invert data
        data = -data;
    end
    
    file_TS_usec = FileStarted_TS(ii_file_start_entry);
 
    % write to Nlx file format

    line = fprintf('%s', temp);
    
    TimestampsSaveflag = 0;
    for cnl = data_cnl_ind
        
        cnl_data = data(cnl+1,:);
        fileName = fullfile(Nlx_OutDir, sprintf('CSC%i.dat', cnl+1));

        channelFid = fopen(fileName, 'a');
        
        fwrite(channelFid, cnl_data, 'int16');
        fclose(channelFid);
        
        if TimestampsSaveflag == 0
            timestampsFilename = fullfile(Nlx_OutDir, 'Timestamps.mat');
            timestampsData = file_TS_usec + (0 : str2double(ADC_period_usec{1}) : str2double(ADC_period_usec{1}) * (numel(cnl_data)-1));
            
            if ii_file_start_entry == 1
                save(timestampsFilename, 'timestampsData');
                
            else
                m = matfile(timestampsFilename, 'Writable', true);
                m.timestampsData = [m.timestampsData timestampsData];
            end
            
            TimestampsSaveflag = 1;
        end
        
%         cnl_data_blocks = vec2mat(cnl_data,512)';
%         num_blocks = size(cnl_data_blocks, 2);
%         blocks_timestamps_usec = file_TS_usec + (0:block_period_time_usec:(block_period_time_usec*(num_blocks-1)));
%         file_name = ['CSC' num2str(cnl) '.ncs'];
%         out_file = fullfile(Nlx_OutDir, file_name);
%         
%         % we do this because Nlx functions have bugs with writing header
%         % when using append... (this way works...)
%         if exist(out_file, 'file')
%             append_flag = 1;
%             Mat2NlxCSC(out_file, append_flag, 1, [], [1 0 0 0 1 0], blocks_timestamps_usec, cnl_data_blocks);
%         else
%             append_flag = 0;
%             Mat2NlxCSC(out_file, append_flag, 1, [], [1 0 0 0 1 1], blocks_timestamps_usec, cnl_data_blocks, header );
%         end
     
    end
    fprintf(repmat('\b', 1, line));
end

end