features_for_header = load('Header_Features.mat');
fet_for_head = strtrim(features_for_header.Header2);

file_search = sprintf('D:\\experiment_data\\*.ntt');
files_to_load = subdir(file_search);
FieldSelection = [1 1 1 1 1];
ExtractionMode = 1;
idx_session = [];

for i_ntt = 1:length(files_to_load)

    current_file = files_to_load(i_ntt).name;
    fprintf('%s\t%.2f\n', current_file, i_ntt*100/length(files_to_load));
    [Timestamps_sessions, ScNumbers, CellNumbersSpikeSorting, Features, Samples, NlxHeader] =...
        Nlx2MatSpike( current_file, FieldSelection, 1, ExtractionMode, idx_session);
    
    CellNumbersSpikeSorting(CellNumbersSpikeSorting ~= 0) = 0;
    
    NlxHeader(end-7:end) = [];
    
    NlxHeader = [NlxHeader; fet_for_head];   
    
    Mat2NlxSpike( current_file, 0, 1, 1,...
        [1 1 1 1 1 1], Timestamps_sessions, ScNumbers, CellNumbersSpikeSorting, Features, Samples, NlxHeader);
end