function p = PRE_get_session_times(p)
%
% prompt the user to get the event timestamps of the session
%

datadir_in = sprintf('%s\\%s\\',p.path_datain, p.data_dir);

[event_list,time_stamps] = get_events_timestamps(datadir_in);

p.event_list = event_list;
p.time_stamps = time_stamps;

nsessions = length(p.S);

% for each session - find the corresponding timestamps of start and end

for nses = 1:nsessions
    s = p.S(nses);
    start_time = time_stamps(s.events(1))+s.time_offsets(1)*1e6;
    end_time = time_stamps(s.events(2))+s.time_offsets(2)*1e6;
    p.S(nses).start_time = start_time;
    p.S(nses).end_time = end_time;
    if s.time_offsets(1) == 0
        start_event = event_list{s.events(1)};
    elseif s.time_offsets(1) > 0
        start_event = sprintf('%s+%d',event_list{s.events(1)},...
            s.time_offsets(1));
    else
        start_event = sprintf('%s-%d',event_list{s.events(1)},...
            -s.time_offsets(1));
    end
    if s.time_offsets(2) == 0
        end_event = event_list{s.events(2)};
    elseif s.time_offsets(2) > 0
        end_event = sprintf('%s+%d',event_list{s.events(2)},...
            s.time_offsets(2));
    else
        end_event = sprintf('%s-%d',event_list{s.events(2)},...
            -s.time_offsets(2));
    end
    p.S(nses).start_event = start_event;
    p.S(nses).end_event = end_event;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EventStrings,Event_TimestampVariable] = get_events_timestamps(datadir_in)

file_name = sprintf('%s\\Events.nev',datadir_in);

%Next we want to extract the Arena Marking Events as well
Behav_events_Extract_Fields = [1 0 0 0 1] ; %We want to extract fields for Timestamps and Event Strings Only.
Extract_Header = 0 ; % Extract the Header as well
Extraction_Mode = 1 ; % Extract All
ExtractionModeArray = [] ; % Extract All

[Event_TimestampVariable, EventStrings] = ...
    Nlx2MatEV( file_name, ...
    Behav_events_Extract_Fields, Extract_Header,Extraction_Mode);

disp('');




