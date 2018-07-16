function p = NLG_PRE_get_session_times_Nlg_and_Nlx(p)


%% In Nlg format - for each session - find the corresponding timestamps
% of its start and end, and update the structure. 
% Both Nlx and Nlg timestamps are in microseconds(!)

%----------------------------------------------------------------------------------------
for nses = 1:length(p.S)
    s = p.S(nses);
    
        start_time = p.Nlg_EventTimestamps(s.events(1))+s.time_offsets_in_seconds(1)*1e6;
        end_time = p.Nlg_EventTimestamps(s.events(2))+s.time_offsets_in_seconds(2)*1e6;
        p.S(nses).start_time = start_time;
        p.S(nses).end_time = end_time;
        
        if s.time_offsets_in_seconds(1) == 0
            start_event = p.Nlg_EventStrings{s.events(1)};
        elseif s.time_offsets_in_seconds(1) > 0
            start_event = sprintf('%s+%d',p.Nlg_EventStrings{s.events(1)},...
                s.time_offsets_in_seconds(1));
        else
            start_event = sprintf('%s-%d',p.Nlg_EventStrings{s.events(1)},...
                -s.time_offsets_in_seconds(1));
        end
        if s.time_offsets_in_seconds(2) == 0
            end_event = p.Nlg_EventStrings{s.events(2)};
        elseif s.time_offsets_in_seconds(2) > 0
            end_event = sprintf('%s+%d',p.Nlg_EventStrings{s.events(2)},...
                s.time_offsets_in_seconds(2));
        else
            end_event = sprintf('%s-%d',p.Nlg_EventStrings{s.events(2)},...
                -s.time_offsets_in_seconds(2));
        end
        
        
        p.S(nses).start_event = start_event;
        p.S(nses).end_event = end_event;
            
        % if we want to  use manually entered timestamp (in microseconds) (good for cases when we are missing the exact event)

        if s.exact_time_in_microseconds(1) ~=0
            p.S(nses).start_time = s.exact_time_in_microseconds(1);
            p.S(nses).start_event=sprintf('User defined %d (microseconds)',p.S(nses).start_time);
        end
        if s.exact_time_in_microseconds(2) ~=0
            p.S(nses).end_time =  s.exact_time_in_microseconds(2);
            p.S(nses).end_event=sprintf('User defined %d (microseconds)',p.S(nses).end_time);
        end;
end

%% In Nlx format - for its session - find the corresponding timestamps
% of its start and end, and update the structure;
% The Nlx timestamps are in microseconds while Nlg ones are in milliseconds
%----------------------------------------------------------------------------------------


file_name = [p.path_datain p.path_year_bat p.path_day_dir '\Nlx_VT_and_TTL\Events.nev']

%Next we want to extract the Arena Marking Events as well
Behav_events_Extract_Fields = [1 0 0 0 1] ; %We want to extract fields for Timestamps and Event Strings Only.
Extract_Header = 0 ; % Extract the Header as well
Extraction_Mode = 1 ; % Extract All
ExtractionModeArray = [] ; % Extract All
[Event_TimestampVariable, EventStrings] = ...
    Nlx2MatEV( file_name, ...
    Behav_events_Extract_Fields, Extract_Header,Extraction_Mode);

%We assume that we want to extact only one session (containing all the data
%including sleep and LED callibrations because we will use it for finding
%the TTLs within this period)
Nlx_timestamps(1)=Event_TimestampVariable(p.Nlx_events(1))+p.Nlx_time_offsets_in_seconds(1)*1e6;
Nlx_timestamps(2)=Event_TimestampVariable(p.Nlx_events(2))+p.Nlx_time_offsets_in_seconds(2)*1e6;
        
% if we want to  use manually entered timestamp (in microsecondss) (good for cases when we are missing the exact event)

if p.Nlx_exact_time_in_microseconds(1) ~=0
    Nlx_timestamps(1)=Nlx_exact_time_in_microseconds(1);
end
if p.Nlx_exact_time_in_microseconds(2) ~=0
    Nlx_timestamps(2)=Nlx_exact_time_in_microseconds(2);
end;
%updating the structure
p.Nlx_timestamps=Nlx_timestamps;


