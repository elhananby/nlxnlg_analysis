function [Samples Timestamps CSC_Sampling_Rate_Hz] = Nlg_CutSamples_CSC(filename, tstart, tend)

%---------------------------
% Arseny Finkelstein (adapted from Michael Yartsev telemetry code) 03/04/2014
%--------------------------

% Cut out data from CSC file according to given timestamps
% Just a wrapper for NLX - Nlx2MatCSC.m
% Input - 
%   filename - full path
%   tstart-tend - in usec
% Output - 
%   Samples - (in uv - assuming that this is what the Neurologger is generating)
%   Timestamps - usec
%

%%  Extract CSCs

%Input variables definition
%------------------------------
filename_CSC_EEG_in = filename;

FieldSelection = [1 1 1 1 1] ; % Will read all the variables from the file

ExtractHeader = 1; %We want to extract the header

ExtractMode = 4; %We want to extract all the records according to given timestamsps
ExtractionModeArray(1) = tstart;
ExtractionModeArray(2) = tend;

[Timestamps, ChanNum, SampleFrequency, NumValidSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, 1, 1, [] ) ;

Timestamps = Timestamps(Timestamps >= tstart & Timestamps <= tend);
Samples = Samples(:, Timestamps >= tstart & Timestamps <= tend);

% for i = 1:length(NlxHeader)
%     if findstr('ADBitVolts',NlxHeader{i})
%         conversion_factor_bits_to_uvolts = str2num(NlxHeader{i}(12:24)) * 10^6
%     end
% end

%% Arseny - add manually (in the correct way - this should be read from the HEADER)
conversion_factor_bits_to_uvolts = 1;
CSC_Sampling_Rate_Hz = regexpi(NlxHeader(contains(NlxHeader, '-SamplingFrequency ')), '(\d+.\d)', 'match');
CSC_Sampling_Rate_Hz = str2double(CSC_Sampling_Rate_Hz{1});

CSC_SamplePeriod_microsec = round( 10^6/CSC_Sampling_Rate_Hz );
Samples_reshaped = zeros(1, numel(Samples) );
Timestamps_filledIn = zeros(1, numel(Samples) );
for ii_DataBlock = 1:size(Samples,2) % Loop over the 512-point blocks of data
    idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
    Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
    Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
end
Samples = Samples_reshaped - mean(Samples_reshaped); % The data, with mean removed
Timestamps = Timestamps_filledIn; % Timestamps in microsec
Samples_corrected_to_uv = conversion_factor_bits_to_uvolts*Samples;
