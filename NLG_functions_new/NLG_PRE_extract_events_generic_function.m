function [Event_TimestampVariable, EventStrings] =	NLG_PRE_extract_events_generic_function(file_name)

% Read the events file to extract the times of the different sessions (Sleep and behavior). 
% The data is not saved from this script but inserted manually into the other scripts which will follow (see below)

% Neuralynx Converters directory must be defined in the matlab path


%------------------------------------------------------------
   
Behav_events_Extract_Fields = [1 0 0 0 1];	% We want to extract fields for Timestamps and Event Strings Only.
Extract_Header				= 0;			% Extract the Header as well
Extraction_Mode				= 1;			% Extract All
ExtractionModeArray			= [];			% Extract All

[Event_TimestampVariable, EventStrings] =	...
    Nlx2MatEV( file_name,					...
    Behav_events_Extract_Fields, Extract_Header,Extraction_Mode);

Event_TimestampVariable=Event_TimestampVariable';
