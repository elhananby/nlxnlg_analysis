function [Samples, Timestamps, CSC_Sampling_Rate_Hz] = Nlg_CutSamples_CSC(filename, tstart, tend)

%Input variables definition
%------------------------------
filename_CSC_EEG_in = filename;

FieldSelection = [1 0 0 0 1] ; % Will read all the variables from the file

ExtractHeader = 1; %We want to extract the header

ExtractMode = 1; %We want to extract all the records according to given timestamsps
ExtractionModeArray(1) = tstart;
ExtractionModeArray(2) = tend;

[Timestamps, Samples, NlxHeader] = ...
    Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, ExtractHeader, ExtractMode, ExtractionModeArray ) ;
timestampsFull = linspace(Timestamps(1), Timestamps(end)+mean(diff(Timestamps)), numel(Samples));
samplesFull = Samples(:);

keepIdx = timestampsFull >= tstart & timestampsFull <= tend;
Samples = samplesFull(keepIdx) - mean(samplesFull(keepIdx));
Timestamps = timestampsFull(keepIdx);
CSC_Sampling_Rate_Hz = str2double(cell2mat(regexp(NlxHeader{8}, '\d*\.\d', 'match')));

end
% for i = 1:length(NlxHeader)
%     if findstr('ADBitVolts',NlxHeader{i})
%         conversion_factor_bits_to_uvolts = str2num(NlxHeader{i}(12:24)) * 10^6
%     end
% end

%% Arseny - add manually (in the correct way - this should be read from the HEADER)
% conversion_factor_bits_to_uvolts=1;


% CSC_SamplePeriod_microsec = round( 10^6/CSC_Sampling_Rate_Hz );
% Samples_reshaped = zeros(1, prod(size(Samples)) );
% Timestamps_filledIn = zeros(1, prod(size(Samples)) );
% for ii_DataBlock = 1:size(Samples,2) % Loop over the 512-point blocks of data
%     idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
%     Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
%     Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
% end
% Samples = Samples_reshaped - mean(Samples_reshaped); % The data, with Mean Removed
% Timestamps = Timestamps_filledIn; % Timestamps in microsec
% Samples_corrected_to_uv = conversion_factor_bits_to_uvolts*Samples;
