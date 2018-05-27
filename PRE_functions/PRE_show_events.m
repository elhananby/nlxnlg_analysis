function PRE_show_events(p, type)
%
% plot out the event timestamps of the session
%
if isequal(type, 'calibration')
    datadir_in = fullfile(p.path_datain, p.calibration_dir);
else
    datadir_in = fullfile(p.path_datain, p.data_dir);
end

[event_list,time_stamps] = get_events_timestamps(datadir_in, type);

fprintf('%s Event list:\n===============\nNum\tString\t\t\t\tTime (min)\n', type);


for i = 1:length(event_list)
    if any(strcmpi(event_list{i}, {'Starting Recording', 'Stopping Recording'}))
        fprintf('%i)\t%s\t%.2f\n', i, event_list{i}, ((time_stamps(i)-time_stamps(1))*1e-6)/60);
    end
end

fprintf('Total Recording Time = %.2f minutes\n', ((time_stamps(end)-time_stamps(1))*1e-6)/60);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EventStrings,Event_TimestampVariable] = get_events_timestamps(datadir_in, eventfile_type)

if isequal(eventfile_type, 'calibration') || isequal(eventfile_type, 'recording')
    file_name = fullfile(datadir_in, 'Events.nev');
    
elseif eventfile_type == 'throw_away'
    file_name = fullfile(datadir_in, 'times_to_throw_away.nev' );
    
end

%Next we want to extract the Arena Marking Events as well
Behav_events_Extract_Fields = [1 0 0 0 1] ; %We want to extract fields for Timestamps and Event Strings Only.
Extract_Header = 0 ; % Extract the Header as well
Extraction_Mode = 1 ; % Extract All
ExtractionModeArray = [] ; % Extract All

[Event_TimestampVariable, EventStrings] = ...
    Nlx2MatEV( file_name, ...
    Behav_events_Extract_Fields, Extract_Header,Extraction_Mode);

end
