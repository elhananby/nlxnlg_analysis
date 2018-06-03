function p = NLG_PRE_process_data(p)
global useGPU;
useGPU = 0;

% 1  convert to CSC
p = NLG_PRE_Nlg2Nlx(p);

% 2 Extract Nlg CSC data, divide into 1-min chunks, and filter it
NLG_filter_CSC(p); % runs on all channels including the 'bad' ones

% 3 Clean CSC Artifacts epochs
NLG_clean_artifacts(p); % runs on all channels including the 'bad' ones

% 4 Detect spikes from CSC
NLG_detect_spikes(p);

% 5 Sync Nlx2Nlg by TTLs
p = NLG_PRE_sync_Nlx2Nlg_by_TTL_2_0(p);
end