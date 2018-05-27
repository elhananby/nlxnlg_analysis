function [ Nlg_EventStrings, Nlg_EventTimestamps, TTL_timestamps_Nlg] = NLG_PRE_read_excel_events_and_TTL(file_name_Event_Nlg)



% current_file_to_read = [file_name_Event_Nlg 'EVENTS.nev'];
% current_file_to_read = [file_name_Event_Nlg 'EVENTS__Mode change.nev'];
% current_file_to_read = [file_name_Event_Nlg 'EVENTS__Parameter change.nev'];
% current_file_to_read = [file_name_Event_Nlg 'EVENTS__PC-generated comment.nev'];
% current_file_to_read = [file_name_Event_Nlg 'EVENTS__Recording parameters.nev'];
% current_file_to_read = [file_name_Event_Nlg 'EVENTS__Clocks synchronized.nev'];

%Extract TTL timestamps
current_file_to_read = [file_name_Event_Nlg 'EVENTS__Digital out.nev'];
[Event_TimestampVariable, EventStrings] = NLG_PRE_extract_events_generic_function (current_file_to_read);
TTL_timestamps_Nlg=Event_TimestampVariable;

%Extract event names
current_file_to_read = [file_name_Event_Nlg 'EVENTS__Free text.nev'];

if exist(current_file_to_read)
    [Event_TimestampVariable, EventStrings] = NLG_PRE_extract_events_generic_function (current_file_to_read);
    Nlg_EventStrings=EventStrings;
    for ii_event=1:1:size(Nlg_EventStrings,1)
        Nlg_EventStrings{ii_event,1}=[ num2str(ii_event) ') ' Nlg_EventStrings{ii_event,1} ];
    end
    
else
    current_file_to_read = [file_name_Event_Nlg 'EVENTS.nev'];
    [Event_TimestampVariable, EventStrings] = NLG_PRE_extract_events_generic_function (current_file_to_read);
    Nlg_EventStrings=EventStrings;
    for ii_event=1:1:size(Nlg_EventStrings,1)
        Nlg_EventStrings{ii_event,1}=[ num2str(ii_event) ') ' Nlg_EventStrings{ii_event,1} ];
    end
end

%Extract event timestamps
Nlg_EventTimestamps=Event_TimestampVariable;

