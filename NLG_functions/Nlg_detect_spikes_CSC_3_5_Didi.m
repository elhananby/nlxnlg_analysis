
function Nlg_detect_spikes_CSC_3_5(expname,forcerecalc)
%
% function Nlg_detect_spikes_CSC_2_1(expname,forcerecalc)
%
% Didi Omer,
% November 2015
%
%
% todo:
%
%

tic1 = tic;
disp('start: Nlg_detect_spikes_CSC_3_5()');
eval(expname);

folder_name = [param.path.analysis, 'preprocessing_spikes3',filesep];
if exist(folder_name,'dir') && ~forcerecalc
    disp(['Detect-spikes was already done',expname]);
    return;
end
thresholds = zeros(4,4); 
tts = reshape([1:16],4,4);
for ii_tetrode=1:length(param.tetrodes.use_for_sorting)
    
    tic2 = tic;
    filename = cell(1,4);
    Spiketimeidx = cell(1,4);
    csc = cell(4,1);
    path = [param.path.CSC_spikes,'spikes_bat',param.bat,'_',param.day,'_TT'];
    tetrodes_use_for_sorting = param.tetrodes.use_for_sorting;
    x_sep_spike_thres= param.spikes.x_sep_spike_thres;
    thresh = zeros(1,4);
   
    % detect spikes
    parfor ch=1:4
        filename= [path,num2str(tts(ch,tetrodes_use_for_sorting(ii_tetrode)))];
        csc{ch,1} = load(filename);
        %thresh = std(csc{ch,1}.Spikes.data.Samples).*4 + (mean(csc{ch,1}.Spikes.data.Samples)); % changed following a disucssion with Nachum and Tamir (29.12.2016)
        thresh(ch) = 2.8* median(abs(csc{ch,1}.Spikes.data.Samples(:))/0.6745); % thresh by Quiroga 2004 (Didi Jan 2017)
        [spikeidx,spikepk]  = peakseek(csc{ch,1}.Spikes.data.Samples,x_sep_spike_thres,thresh(ch));
        [spikeidx,ii] = sort(spikeidx);
        spikepk = spikepk(ii);
        spikeidx = spikeidx((spikeidx>7 & spikeidx<=length(csc{ch,1}.Spikes.data.Samples)-24));
        spiketime = csc{ch,1}.Spikes.data.Timestamps(spikeidx);
        spiketimeidx{ch,1} = [spiketime(:)';spikeidx(:)'];
       
    end
    threshold(:,ii_tetrode)= thresh(:); 
    SamplingFreq = csc{ii_tetrode}.Spikes.params.CSC_SamplingFreq;
    disp(['1) detecting spikes ',expname,' on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic2)/60),' min']);
    
    %% coincedent detection (matrixwise implementation) and removal of spikes which fall within a refractory period
    tic3 = tic;
    startidx = min([spiketimeidx{1}(2,:)' ;spiketimeidx{2}(2,:)';spiketimeidx{3}(2,:)';spiketimeidx{4}(2,:)']);
    endidx  =  max([spiketimeidx{1}(2,:)' ;spiketimeidx{2}(2,:)';spiketimeidx{3}(2,:)';spiketimeidx{4}(2,:)']);
    len = startidx - endidx +1;
    M = zeros(4,len);
    for i=1:4
        M(i,spiketimeidx{i}(2,:))=csc{i,1}.Spikes.data.Samples(spiketimeidx{i}(2,:));
    end
    M = max(M);
    
    spikeidx = peakseek(M,22,min(thresh(:)));
    spiketime = csc{1,1}.Spikes.data.Timestamps(spikeidx);
    clear M;
    clear xx;
    clear idx;
    disp(['2) remove redundent spike events',expname,' on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic3)/60),' min']);
    
    %% calculate spikewave
    tic4 = tic;
    spikewave = zeros(32,size(spiketime,2),4);
    mat= repmat([-7:24]',1,length(spikeidx));
    spike_idx = repmat(spikeidx,[32,1]);
    spike_idx = spike_idx+mat;
    parfor ch=1:4
        spikewave(:,:,ch) = csc{ch,1}.Spikes.data.Samples(spike_idx);
    end
    disp(['3) calculating waveforms TT',num2str(ii_tetrode),' (',num2str(size(spiketime,2)),' spikes)',' elepsed:',num2str(toc(tic4)/60),' min']);
    clear csc;
    
    %% compare to acceptable spike shapes
    tic6 = tic;
    eval(['load ', param.path.spike_shapes_lib]); % Load the library of acceptable spike shapes;
    [~,iimax ] = max(squeeze(max(spikewave)),[],2);
    waves = zeros(size(spikewave,1),size(spikewave,2));
    for i=1:size(waves,2)
        waves(:,i) = spikewave(:,i,iimax(i));
    end
    library_of_acceptable_spike_shapes = library_of_acceptable_spike_shapes';
    c = zeros(size(library_of_acceptable_spike_shapes,2),size(waves,2));
    parfor ii=1:size(library_of_acceptable_spike_shapes,2)
        mat = repmat(library_of_acceptable_spike_shapes(:,ii),1,size(waves,2));
        c1 = corrcoeff(waves(2:end-1,:),mat(2:end-1,:));
        c2 = corrcoeff(waves(1:end-2,:),mat(2:end-1,:));
        c3 = corrcoeff(waves(3:end,:),mat(2:end-1,:));
        c(ii,:) = max([c1;c2;c3]);
    end
    c = max(c);
    acceptedidx = find(c>=0.8);
    Timestamps_accepted_spikes_TT{ii_tetrode}=spiketime(acceptedidx);
    spikes_TT{ii_tetrode}=spikewave(:,acceptedidx,:);
    disp(['4) Testing spikewaves against spike-shape lib. TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic6)/60),' min' ]);
    clear c;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Coincidence-Detection across Tetrodes to eliminate artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic7 = tic; 
if length(param.tetrodes.use_for_sorting)>1
    idx_coincidence_vec = cell(length(param.tetrodes.use_for_sorting),1);
    
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}), % Loop over the spikes of the FIRST tetrode
        temp_stack = cell(length(param.tetrodes.use_for_sorting),1);
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        temp_stack{1,1} = ii_spikes;
        for TT=2:length(param.tetrodes.use_for_sorting)
            temp_stack{TT,1} = find( abs( Timestamps_accepted_spikes_TT{TT} - t_spike ) <= param.spikes.coincidence_window ); % THE COINCIDENCE DETECTION
        end
        ii = ~cellfun(@isempty,temp_stack);
        if sum(ii)>=length(param.tetrodes.use_for_sorting)
            for jj=1:length(param.tetrodes.use_for_sorting)
                idx_coincidence_vec{jj,1} = cat(2,idx_coincidence_vec{jj,1}, temp_stack{jj,1});
            end
        end
    end
end
   

for i=1:length(param.tetrodes.use_for_sorting)
    idx_coincidence_vec{i} = unique(idx_coincidence_vec{i});
    if ~isempty(idx_coincidence_vec{i})
        Timestamps_accepted_spikes_TT{i}(idx_coincidence_vec{i})=[];
        spikes_TT{i}(:,idx_coincidence_vec{i},:)=[];
    end
end
disp(['5) cross TT coincident spike removal...',' elepsed:',num2str(toc(tic7)/60),' min']);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
%%%%%%%%%%%%%%%%%%%%%%%%%

tic8 = tic;
for ii_Tetrode=1: length(param.tetrodes.use_for_sorting)   
          
    base_output_dir						= [param.path.analysis, 'preprocessing_spikes4',filesep];
    if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end
    filename_out = [base_output_dir, 'spikes_',param.DirRaw,'_TT', num2str(param.tetrodes.use_for_sorting(ii_Tetrode))];
    
    % constract header for NTT file
    load(param.path.NTT_header);
    l = ['-ThreshVal ',num2str(threshold(1,ii_Tetrode)),' ',...
                    num2str(threshold(2,ii_Tetrode)),' ',...
                    num2str(threshold(3,ii_Tetrode)),' ',...
                    num2str(threshold(4,ii_Tetrode))]; 
    Header{end+1} = l;
    Header = MakeNlxHeader(Header);
    %Header = Header'; 
    
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{ii_Tetrode};
    spikes=spikes_TT{ii_Tetrode};
    spikes = permute(spikes,[1,3,2]);
    
    %only export timestamps and data points - Full session
    FieldSelection = [1 0 0 0 1 1];
    
    if ispc
        Mat2NlxSpike([filename_out,'.NTT'],0,1,[],FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    elseif isunix
        if exist([filename_out,'.NTT'],'file'), delete([filename_out,'.NTT']),end
        Mat2NlxTT([filename_out,'.NTT'],0,1,1,length(Timestamps_accepted_spikes(:)),FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    end
end
disp(['6) Saving spikes to disk...',' elepsed:',num2str(toc(tic8)/60),' min']);
disp(['8) Total detection on all TTs ',expname,' elepsed:',num2str(toc(tic1)/60),' min']);

end



