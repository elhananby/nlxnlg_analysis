function p = NLG_PRE_process_data(p)
% 1  convert to CSC
NLG_PRE_Nlg2Nlx(p);
% NLG_Nlg2Nlx(p);

% 2 Extract Nlg CSC data, divide into 1-min chunks, and filter it
Nlg_Flight_extract_and_filter_CSC_parfor(p); % runs on all channels including the 'bad' ones
% NLG_Filter_Channel(p);

% 3 Clean CSC Artifacts epochs
% Nlg_clean_artifacts_CSC_parfor(p); % runs on all channels including the 'bad' ones

% 4 Detect spikes from CSC
% Nlg_Flight_detect_spikes_CSC_parfor(p); %CSC artifacts are removed from the neural data based on 'good' channels only
NLG_detect_spikes(p);

% 5 Sync Nlx2Nlg by TTLs
p = NLG_PRE_sync_Nlx2Nlg_by_TTL_2_0(p);
end