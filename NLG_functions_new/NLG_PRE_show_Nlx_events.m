function  [EventStrings Event_TimestampVariable] = NLG_PRE_show_Nlx_events(p)

%---------------------------
% Arseny Finkelstein 18/04/2014

% display the event timestamps of the session
%--------------------------


datadir_in = [p.path_datain p.path_year_bat p.path_day_dir];

file_name = sprintf('%s\\Nlx_VT_and_TTL\\Events.nev',datadir_in);

%Next we want to extract the Arena Marking Events as well    
Behav_events_Extract_Fields = [1 0 0 0 1] ; %We want to extract fields for Timestamps and Event Strings Only.
Extract_Header = 0 ; % Extract the Header as well
Extraction_Mode = 1 ; % Extract All
ExtractionModeArray = [] ; % Extract All

[Event_TimestampVariable, EventStrings] = ...
    Nlx2MatEV( file_name, ...
    Behav_events_Extract_Fields, Extract_Header,Extraction_Mode);


for i = 1:length(EventStrings)
    disp(strcat(num2str(i),') ',EventStrings{i}));
end
disp(' ');

