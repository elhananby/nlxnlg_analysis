function  NLG_DBG_find_TTL_delay_by_LED (file_name_Nlx_TTL,folder_name, p)
% =======================================================================

%---------------------------
% Arseny Finkelstein 18/03/2014

% Finding the delay introduced between TTL_pulse and TTL-out event in the Neurologger. We do it by synchronizing with LEDs on/off times (in case we don't have TTLs)

%--------------------------

filename_data_to_save=[folder_name '\Nlx_VT_and_TTL\TTLs_and_Nlx2Nlg_clock_fits.mat'];


%% 1)  Get Nlg and Nlx sync times 

%Getting Nlg sync times
disp('Showing NLG events for LED sync')
p.Nlg_EventStrings
Nlg_LED_turn_on_events=[2;5;8;11;14;17;20;23;26;29;32;35;38;41;44;46;49;52;55;58];
Nlg_LED_turn_off_events=Nlg_LED_turn_on_events+1;
Nlg_LED_turn_on_times=cell2mat(p.Nlg_EventTimestamps(Nlg_LED_turn_on_events))'*1000; %in microsec
Nlg_LED_turn_off_times=cell2mat(p.Nlg_EventTimestamps(Nlg_LED_turn_off_events))'*1000;%in microsec

%Getting Nlx sync times
disp('Showing NLG events for LED sync')
p.Nlx_EventStrings
Nlx_LED_turn_on_events=[2:2:40]';
Nlx_LED_turn_off_events=Nlx_LED_turn_on_events+1;
Nlx_LED_turn_on_times=p.Nlx_EventTimestamps(Nlx_LED_turn_on_events); %in microsec
Nlx_LED_turn_off_times=p.Nlx_EventTimestamps(Nlx_LED_turn_off_events);%in microsec

% Extract the video files from both cameras for synchronization:
% --------------------------------------------------
VT_directory=[folder_name '\Nlx_VT_and_TTL\'];

for ii_sync=1:1:length(Nlx_LED_turn_on_times)
    
    Extract_Fields = [1 1 1 0 1 0] ; % For each frame extract all variables except ExtractedAngle & Points
    Extract_Header = 1 ; % Extract the Header as well
    Extraction_Mode = 4 ; % Extract a subset of video frames, specified by an Timestamp Range
    start_Nlx_sync = Nlx_LED_turn_on_times (ii_sync);
    end_Nlx_sync = Nlx_LED_turn_off_times (ii_sync);
    ExtractionModeArray = [start_Nlx_sync end_Nlx_sync]; % Timestamp Range to Extract (it is defined above)
    
    VT1_file = fullfile(VT_directory,'VT1.nvt');
    [Timestamps_v1, ExtractedX, ExtractedY, Targets_VT1, NlxHeader] = ...
        Nlx2MatVT( VT1_file, Extract_Fields, Extract_Header,Extraction_Mode,ExtractionModeArray) ; % Extract data - camera 1
    
    Nlx_LED_on_Timestamps_v1 (ii_sync)=Timestamps_v1(find(ExtractedX,1,'first'));
%     Nlx_LED_off_Timestamps_v1(ii_sync)=Timestamps_v1(find(ExtractedX,1,'last')+1);
    
    VT2_file = fullfile(VT_directory,'VT2.nvt');
    [Timestamps_v2, ExtractedX, ExtractedY, Targets_VT2, NlxHeader] = ...
        Nlx2MatVT(VT2_file, Extract_Fields, Extract_Header, Extraction_Mode,ExtractionModeArray ) ; % Extract data - camera 2
    
    Nlx_LED_on_Timestamps_v2(ii_sync)=Timestamps_v2(find(ExtractedX,1,'first'));
%     Nlx_LED_off_Timestamps_v2(ii_sync)=Timestamps_v2(find(ExtractedX,1,'last')+1);
    
end;

Nlx_LED_on_Timestamps_min_v1_v2=min([Nlx_LED_on_Timestamps_v1;Nlx_LED_on_Timestamps_v2])
% Nlx_LED_off_Timestamps_min_v1_v2=min([Nlx_LED_off_Timestamps_v1;Nlx_LED_off_Timestamps_v2])

[val idx_LED_tmestamp]=sort((Nlg_LED_turn_on_times-Nlx_LED_on_Timestamps_min_v1_v2),'descend')
best=2;
LED_timestamps_Nlg = [Nlg_LED_turn_on_times(sort(idx_LED_tmestamp(1:best),'ascend'))];
LED_timestamps_Nlx = [Nlx_LED_on_Timestamps_min_v1_v2(sort(idx_LED_tmestamp(1:best),'ascend'))];

%% 2) Calculating the  fit (for conversion of Nlx timestamps into Nlg timestamps) and residuals 
polyfit_Nlg2Nlx_microsec = polyfit(LED_timestamps_Nlg, LED_timestamps_Nlx, 1);
polyfit_Nlx2Nlg_microsec = polyfit(LED_timestamps_Nlx, LED_timestamps_Nlg, 1)
clock_lag_Nlg2Nlx_micro_sec = polyval(polyfit_Nlg2Nlx_microsec,0); %if positive - it means that Nlx is ahead of neurologger

% Computing the residuals on final TTL - i.e. based on Good TTL
LED_residuals_times =LED_timestamps_Nlg-  polyval(polyfit_Nlx2Nlg_microsec,LED_timestamps_Nlx); 


% Plotting the fit
figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[21 2 25 15]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[-20.9 -0.6 0 0]);

subplot(1,2,1)
hold on;
plot(LED_timestamps_Nlg,LED_timestamps_Nlx,'ob');
plot ([LED_timestamps_Nlg(1), LED_timestamps_Nlg(end)], [polyval(polyfit_Nlg2Nlx_microsec,LED_timestamps_Nlg(1)) , polyval(polyfit_Nlg2Nlx_microsec,LED_timestamps_Nlg(end))], '-k');
xlabel('Nlg timestamp (micro s)');
ylabel('Nlx timestamp (micro s)');
title(sprintf(' Synchonization: Fit of Nlx and Nlg clocks based on LED times fit \n'));

%Plotting the residuals 
subplot(1,2,2)
hold on;
plot(1:length(LED_residuals_times),LED_residuals_times/1000,'ob');
xlabel('# of LED turn on/off events');
ylabel('Time difference (ms)');
title(sprintf('Residuals of the LEDs based fit \n'));


% Saving the Figure
filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\LED_based_Nlg_Nlx_clocks_sync_and_Residuals'];
%save the fig as .tiff
print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
% %save the fig as .fig
saveas(gca,filename_fig_to_save,'fig');


%% 3) Saving the fit

% Saving the fit 
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlg2Nlx_microsec = polyfit_Nlg2Nlx_microsec;
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlx2Nlg_microsec = polyfit_Nlx2Nlg_microsec;

save(filename_data_to_save,'-struct','TTLs_and_Nlx2Nlg_clock_fits');

%% Parameters
run_in_chunks_to_save_memory=0; % 1 - runs cross-corr in chunks to avoid memory overload or 0 - run in one piece (for stronger computers)
heaviside_template_window=50; %in number of samples. Thi is the half window size of the step function (before the step, and the same length after the step).
cross_corr_threshold=5*10^5; %used to identify the peaks of the crosscorrelation, corresponding to the rising phase of the TTL
cross_corr_detection_peak_shift=50; %in units of samples (for heaviside template of length 100(50 -1 and 50 +1);
max_jitter_between_TTL_pairs=80;  %in ms  - this is the maximal acceptablle jitter between the Nlx and Nlg TTLs pairs, that will be used for synchronization
max_residual_TTL_fit_ms=15; % if the residual time between the data and the fit is bigger than this, we will recalculate the fit without this TTL time
filename_data_to_save=[folder_name '\Nlx_VT_and_TTL\TTLs_and_Nlx2Nlg_clock_fits.mat'];
error_sync=0; %intializing. If there will be an error in sync (i.e. no good TTLs pairs) it will change to 1

%% ======== Extract CSC min-max data  ============

% skip CSC Artifact cleaning if we had done this already
if exist ((filename_data_to_save))==2
    disp('======================');
    disp (['Skipping Extracting Nlx TTLs for. Already extracted in ' filename_data_to_save]);
    disp('======================');
    load(filename_data_to_save);
    return;
else
    % =======================================================================
    disp('===========');
    disp(['Extracting Nlx TTLs for ',file_name_Nlx_TTL])
    disp('===========');
end;

%We want to extract fields for Timestamps and Event Strings Only.
FieldSelection(1) = 1;
FieldSelection(2) = 1;
FieldSelection(3) = 1;
FieldSelection(4) = 1;
FieldSelection(5) = 1;

ExtractHeader = 1; %We want to extract the header

ExtractMode = 1; %We want to extract all the records from the file into matlab.
%Extraction Mode Array is blank because we want to extract All, so this mode does not
%require an Array, this parameter should be left blank
ExtractionModeArray=[];

% ExtractMode = 4; %We want to extract all the records within a timestamp range.
% ExtractionModeArray(1) = tstart;
% ExtractionModeArray(2) = tend;
[Timestamps, ChanNum, SampleFrequency, NumValidSamples, Samples, CSC_HeaderVariable] = ...
    Nlx2MatCSC(file_name_Nlx_TTL, FieldSelection, ExtractHeader, ExtractMode, ExtractionModeArray) ;

CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
CSC_Sampling_Rate_Hz = 10^6 / CSC_SamplePeriod_microsec ;
Samples_reshaped = zeros(1, prod(size(Samples)) );
Timestamps_filledIn = zeros(1, prod(size(Samples)) );
for ii_DataBlock = 1:size(Samples,2), % Loop over the 512-point blocks of data
    idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
    Samples_reshaped( idx_data ) = Samples(:,ii_DataBlock)';
    Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
end
Samples_without_DC=Samples_reshaped-mean(Samples_reshaped);

% Extract Nlg TTL times
TTL_timestamps_Nlg = p.TTL_timestamps_Nlg'; %in ms




%% 6) Computing the 'PSTH' for all NLG emitted TTLs synchronized by LED callibration fit

PSTH_window_before_onset= round(0.1*CSC_Sampling_Rate_Hz); %0.1 sec before
PSTH_window_after_onset= round(0.2*CSC_Sampling_Rate_Hz); % %0.5 sec before
TTL_timestamps_Nlx_sync_by_LED= polyval(polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg*1000);
% FOR ALL TTLs
for ii_TTL=1:1:length(TTL_timestamps_Nlx_sync_by_LED)
    [temp TTL_onset_idx_by_LED]=min(abs(Timestamps_filledIn - TTL_timestamps_Nlx_sync_by_LED(ii_TTL)));
    TTL_PSTH(:,ii_TTL) = Samples_reshaped( TTL_onset_idx_by_LED-PSTH_window_before_onset : TTL_onset_idx_by_LED+PSTH_window_after_onset );
end;

% Plotting the 'PSTH' of ALL TTL pulses, aligned to the detection point
figure;
subplot(1,1,1);
hold on;
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),TTL_PSTH,'-k'); %PSTH (all TTL events)
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),mean(TTL_PSTH,2),'-b'); % average TTL shape
plot(1000*[ PSTH_window_before_onset+1 PSTH_window_before_onset+1]./CSC_Sampling_Rate_Hz,[min(min(TTL_PSTH)) max(max(TTL_PSTH))],'-r');
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('All Emitted TTLs aligned to Nlg event time converted to Nlx time by LED on/off sync\n'));


% Saving the Figure
filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\4_TTL_PSTH_aligned_to_peak_detected'];
%save the fig as .tiff
print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
% %save the fig as .fig
saveas(gca,filename_fig_to_save,'fig');












%------------------------------------------
%% 1) Extracting the TTLs peaks by computing crosscorrelation with a step function template
heaviside_template=[-ones(1,heaviside_template_window),ones(1,heaviside_template_window)]; %heaviside step_function imitating the rise time of the TTL pulse

%performing crosscorrelation between the data and the step_function

if run_in_chunks_to_save_memory == 1 % if = 1 - runs cross-corr in chunks to avoid memory overload
    
    chunk_window=10^6; % samples
    chunk_start=1;
    TTL_onset_idx=[];
    for ii_chunk=1:1:ceil(length(Samples_without_DC)/chunk_window)
        chunk_end = chunk_start + chunk_window -1;
        if chunk_end>=length(Samples_without_DC) %if we are at the last chunk
            chunk_window = length(Samples_without_DC) - chunk_start; %make the chunk window shorter
            chunk_end = chunk_start + chunk_window -1;
        end
        
        [c,  lags] = xcorr(Samples_without_DC(chunk_start:chunk_end), heaviside_template);
        c=c(chunk_window:end);
        chunk_start = chunk_start + chunk_window; % update chunk start

        %Identify the peaks of the crosscorrelation, corresponding to the rising
        %phase of the TTL pulse
        c_thresholded=c;
        c_thresholded(find(c_thresholded<cross_corr_threshold))=0;
        [ymax,idx_max,ymin,idx_min] = extrema(c_thresholded); %finding the peaks
        idx_max_sorted=sort(idx_max) +  chunk_window*(ii_chunk-1); % making this idx relative to the complete trace
        idx_max_sorted = idx_max_sorted + cross_corr_detection_peak_shift; %shifting the peaks by this empirical shift so they will be exactly aligned to the rise phase of the TTL;
        TTL_onset_idx= [TTL_onset_idx idx_max_sorted];
        
    end % end looping over chunks
    
    TTL_timestamps_Nlx=Timestamps_filledIn(TTL_onset_idx);
    
elseif run_in_chunks_to_save_memory == 0  % if = 0 - run cross-corr in one piece (for stronger computers)
    
    [c,  lags] = xcorr(Samples_without_DC,heaviside_template);
    c=c(length(Samples_reshaped):end);
    
    %Identify the peaks of the crosscorrelation, corresponding to the rising
    %phase of the TTL pulse
    c_thresholded=c;
    c_thresholded(find(c_thresholded<cross_corr_threshold))=0;
    [ymax,idx_max,ymin,idx_min] = extrema(c_thresholded); %finding the peaks
    
    idx_max_sorted=sort(idx_max);
    
    TTL_onset_idx=idx_max_sorted+cross_corr_detection_peak_shift; %shifting the peaks by this empirical shift so they will be exactly aligned to the rise phase of the TTL
    TTL_timestamps_Nlx=Timestamps_filledIn(TTL_onset_idx);

    % Plotting the crosscorrelations peaks and the resulting peaks on the actual data
    figure;
    % Some WYSIWYG options:
    set(gcf,'DefaultAxesFontSize',7);
    set(gcf,'DefaultAxesFontName','helvetica');
    set(gcf,'PaperUnits','centimeters','PaperPosition',[21 2 28 19.5]);
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[-20.9 -0.6 0 0]);
    hold on
    plot(Timestamps_filledIn,c,'-b');
    plot(Timestamps_filledIn,Samples_without_DC*10,'-r');
    plot(Timestamps_filledIn(idx_max),ymax,'b*'); %crosscorrelations peaks
    plot(TTL_timestamps_Nlx,Samples_without_DC(TTL_onset_idx)*10,'*r'); % resultingCSC peaks
    plot([Timestamps_filledIn(1) Timestamps_filledIn(end)],[cross_corr_threshold cross_corr_threshold],'-k'); % Crosscorrelation detection threshold
    title(sprintf('Cross-correlation.    %d identified Nlx TTLs out of %d Nlg TTLs',length(TTL_timestamps_Nlx),length(TTL_timestamps_Nlg)), 'FontSize', 12);
    legend('Crosscorrelation','CSC with TTLs (Audio 2)','Crosscorrelations peaks', 'CSC peaks','Detection Threshold','Location','SouthEast');
    xlabel('Timestamps in Neuralynx format (microseconds)', 'FontSize', 12,  'VerticalAlignment', 'top');
    xlim([Timestamps_filledIn(1) Timestamps_filledIn(end)]);
    
    % Saving the Figure
    filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\1_TTL_peaks_from_crosscorr'];
    %save the fig as .tiff
    print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
    % %save the fig as .fig
    % saveas(gca,filename_fig_to_save,'fig');
    close all;
    
end;










