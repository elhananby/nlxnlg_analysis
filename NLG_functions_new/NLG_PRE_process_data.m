function p = NLG_PRE_process_data(p)
global useGPU;
useGPU = 0;

% 1 convert to CSC
p = NLG_PRE_Nlg2Nlx(p);

% 2 extract and filter CSC data
NLG_extract_and_filter_CSC(p);

% 3 Clean CSC Artifacts epochs
Nlg_clean_artifacts_CSC_parfor(p); %runs on all channels including the 'bad' ones

% 4 detect spikes
Nlg_Flight_detect_spikes_CSC_parfor(p); %CSC artifacts are removed from the neural data based on 'good' channels only

% 5 sync timestamps

end