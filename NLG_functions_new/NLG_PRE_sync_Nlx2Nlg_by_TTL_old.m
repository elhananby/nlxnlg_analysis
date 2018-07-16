function [TTL_timestamps_Nlx  polyfit_Nlg2Nlx_microsec polyfit_Nlx2Nlg_microsec error_sync]= NLG_PRE_sync_Nlx2Nlg_by_TTL (file_name_Nlx_TTL,folder_name,p)
% =======================================================================

%---------------------------
% Arseny Finkelstein 04/02/2014
%--------------------------

% Extracting TTL onset times fro Nlx, by crosscorrelation-threshold
% detection of the continiously sampled TTL channel ('Audio2')

%% Parameters
run_in_chunks_to_save_memory=1; % 1 - runs cross-corr in chunks to avoid memory overload or 0 - run in one piece (for stronger computers)
% Nlg_event_delay_ms=35; % This is the average delay introduced by the PC (between the emission of TTL pulse and the registration of the TTL-out event). Positve value - means that the TTL pulse was emitted before the TTL-out event was saved in the log-file
Nlg_event_delay_ms=0; % This is the average delay introduced by the PC (between the emission of TTL pulse and the registration of the TTL-out event). Positve value - means that the TTL pulse was emitted before the TTL-out event was saved in the log-file

heaviside_template_window=50; %in number of samples. Thi is the half window size of the step function (before the step, and the same length after the step).
cross_corr_threshold=5*10^5; %used to identify the peaks of the crosscorrelation, corresponding to the rising phase of the TTL
cross_corr_detection_peak_shift=50; %in units of samples (for heaviside template of length 100(50 -1 and 50 +1);
max_jitter_between_TTL_pairs=80;  %in ms  - this is the maximal acceptablle jitter between the Nlx and Nlg TTLs pairs, that will be used for synchronization
max_residual_TTL_fit_ms=15; % if the residual time between the data and the fit is bigger than this, we will recalculate the fit without this TTL time
filename_data_to_save=[folder_name '\Nlx_VT_and_TTL\TTLs_and_Nlx2Nlg_clock_fits.mat'];
filename_error=[folder_name '\Nlx_VT_and_TTL\error_sync.mat'];

%intializing.
error_sync=0;  %If there will be an error in sync (i.e. no good TTLs pairs) it will change to 1
TTL_timestamps_Nlx=[];
polyfit_Nlg2Nlx_microsec=[];
polyfit_Nlx2Nlg_microsec=[];

%% ======== Extract CSC min-max data  ============

% skip CSC Artifact cleaning if we had done this already
if exist ((filename_data_to_save))==2
    disp('======================');
    disp (['Skipping Extracting Nlx TTLs. Already extracted in ' filename_data_to_save]);
    disp('======================');
    load(filename_data_to_save);
    return;
elseif  exist ((filename_error))==2
    disp('===========');
    disp(['ERROR in SYNC of Nlx TTLs for ',file_name_Nlx_TTL])
    disp('===========');
    error_sync=1;
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

% ExtractMode = 1; %We want to extract all the records from the file into matlab.
% %Extraction Mode Array is blank because we want to extract All, so this mode does not
% %require an Array, this parameter should be left blank
% ExtractionModeArray=[];

ExtractMode = 4; %We want to extract all the records within a timestamp range.
ExtractionModeArray(1) = p.Nlx_timestamps (1);
ExtractionModeArray(2) = p.Nlx_timestamps (2);
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
TTL_timestamps_Nlg = (p.TTL_timestamps_Nlg'/1000)-Nlg_event_delay_ms; %in ms

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

%% 2) Detecting missalignments (misses or false positive) between Nlx and Nlg TTLs
%Creating a vector from the first to the last TTL,  with bin size of 1 ms,
% in order to see if we missed some event, and if there is a shift

time_vector_TTL_Nlx = 0*[1:1: ceil( (TTL_timestamps_Nlx(end) - TTL_timestamps_Nlx(1))/(10^3))];
time_vector_TTL_Nlg =  0*[1:1: ceil( (TTL_timestamps_Nlg(end) - TTL_timestamps_Nlg(1))/(10^0))];
% Placing 1 at a timestamp (with bin of 1 ms) where a TTL was sent (Nlg) or detected (Nlx)
time_vector_TTL_Nlx (ceil( (TTL_timestamps_Nlx - TTL_timestamps_Nlx(1))/(10^3)) +1 )=1;
time_vector_TTL_Nlg (ceil( (TTL_timestamps_Nlg - TTL_timestamps_Nlg(1))/(10^0)) +1 )=1;

%performing crosscorrelation between the TTL timestamps in Nlx and Nlg
%formats to see if we had 'misses' or 'false positives'
time_vector_TTL_Nlg_smoothed=filtfilt(hamming(10),1,time_vector_TTL_Nlg); %we do a gaussian smoothing to be robust to potential jitter
time_vector_TTL_Nlx_smoothed=filtfilt(hamming(10),1,time_vector_TTL_Nlx); %we do a gaussian smoothing to be robust to potential jitter
[c_Nlg_Nlx, lags] = xcorr(time_vector_TTL_Nlg_smoothed,time_vector_TTL_Nlx_smoothed);
%finding the delay in ms
[temp shift_Nlg_Nlx_ms ]=max(c_Nlg_Nlx);
shift_Nlg_Nlx_ms= shift_Nlg_Nlx_ms- length(time_vector_TTL_Nlg);

% Plotting crosscorrelation between the TTL timestamps in Nlx and Nlg
% and showing the delay if we had it

figure;
hold on;
plot(c_Nlg_Nlx,'-k'); %cross correlation
plot([ length(time_vector_TTL_Nlg), length(time_vector_TTL_Nlg)],[min(c_Nlg_Nlx) max(c_Nlg_Nlx)],'-.r');
legend('Nlg Nlx crosscorr','0 lag');
text(length(time_vector_TTL_Nlg)/10, max(c_Nlg_Nlx), sprintf('Computed peak shift = %d ms', shift_Nlg_Nlx_ms));
xlabel('Time (ms)');

% Saving the Figure
filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\2_TTL_Nlg_Nlx_crosscorr'];
%save the fig as .tiff
print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
% % %save the fig as .fig
% saveas(gca,filename_fig_to_save,'fig');

%% 3) Aligning the TTLs (in case there were misses) and taking only the TTLs that were present both in Nlg and Nlx
% Updating the Nlx timestamps by the shift (in micro s), which could result
% from missalignment of the TTL trains i.e. - false positive/negative)
TTL_timestamps_Nlx=TTL_timestamps_Nlx+shift_Nlg_Nlx_ms*10^3;

%matching the TTLs according to temporal difference
idx_of_TTL_Nlg=find(time_vector_TTL_Nlg==1);
idx_of_TTL_Nlx=find(time_vector_TTL_Nlx==1);
for ii_TTL = 1:1:length(idx_of_TTL_Nlx)
    [abs_temp_diff_to_closest_Nlg_TTL(ii_TTL) idx_closest_Nlg_TTL_for_each_Nlx_TTL(ii_TTL) ]=min(abs((idx_of_TTL_Nlx(ii_TTL)+shift_Nlg_Nlx_ms) - idx_of_TTL_Nlg));
end;
[temp_diff_sorted  idx_temp_diff_sorted ]=sort(abs_temp_diff_to_closest_Nlg_TTL,'ascend');

%taking the best matching pairs
number_of_pairs=min(length(idx_of_TTL_Nlg),length(idx_of_TTL_Nlx));
idx_matching_Nlx_TTL=idx_temp_diff_sorted(1:number_of_pairs);

% Exluding those pairs with temporal jitter of more than
% max_jitter_between_TTL_pairs in ms  -  the maximal acceptablle jitter between the Nlx and Nlg TTLs pairs, that will be used for synchronization
idx_matching_Nlx_TTL =   sort(idx_matching_Nlx_TTL(temp_diff_sorted(1:number_of_pairs)<max_jitter_between_TTL_pairs));
idx_matching_Nlg_TTL =   idx_closest_Nlg_TTL_for_each_Nlx_TTL (idx_matching_Nlx_TTL);

%% 4) Computing the initial fit (for conversion of Nlx timestamps into Nlg timestamps) and residuals on matching TTLs

initial_polyfit_Nlg2Nlx_microsec = polyfit(TTL_timestamps_Nlg(idx_matching_Nlg_TTL)*1000, TTL_timestamps_Nlx(idx_matching_Nlx_TTL), 1);
initial_polyfit_Nlx2Nlg_microsec = polyfit(TTL_timestamps_Nlx(idx_matching_Nlx_TTL),TTL_timestamps_Nlg(idx_matching_Nlg_TTL)*1000,  1);

% Computing the residuals on matching TTLs ( i.e. how much each actual measured TTL was different from the fit)
matching_TTL_residuals_times =TTL_timestamps_Nlg(idx_matching_Nlg_TTL)-  polyval(initial_polyfit_Nlx2Nlg_microsec,TTL_timestamps_Nlx(idx_matching_Nlx_TTL))/1000; %if positive - it means that Nlx is ahead of neurologger

% Identifying 'Good' and 'Bad' TTL pairs based on the residuals of the inital fit
%good
idx_good_TTL_Nlg=idx_matching_Nlg_TTL(find(abs(matching_TTL_residuals_times)<max_residual_TTL_fit_ms));
idx_good_TTL_Nlx=idx_matching_Nlx_TTL(find(abs(matching_TTL_residuals_times)<max_residual_TTL_fit_ms));
%bad
idx_bad_TTL_Nlg=idx_matching_Nlg_TTL(find(abs(matching_TTL_residuals_times)>max_residual_TTL_fit_ms));
idx_bad_TTL_Nlx=idx_matching_Nlx_TTL(find(abs(matching_TTL_residuals_times)>max_residual_TTL_fit_ms));


%% 5) Recalculating the final fit (for conversion of Nlx timestamps into Nlg timestamps) and residuals on Good TTL, i.e. - those TTL whose residuals didn't deviate significantly from the intial fit
polyfit_Nlg2Nlx_microsec = polyfit(TTL_timestamps_Nlg(idx_good_TTL_Nlg)*1000, TTL_timestamps_Nlx(idx_good_TTL_Nlx), 1);
polyfit_Nlx2Nlg_microsec = polyfit(TTL_timestamps_Nlx(idx_good_TTL_Nlx),TTL_timestamps_Nlg(idx_good_TTL_Nlg)*1000,  1);
clock_lag_Nlg2Nlx_micro_sec = polyval(polyfit_Nlg2Nlx_microsec,0); %if positive - it means that Nlx is ahead of neurologger

% Computing the residuals on final TTL - i.e. based on Good TTL
good_TTL_residuals_times =TTL_timestamps_Nlg(idx_good_TTL_Nlg)-  polyval(polyfit_Nlx2Nlg_microsec,TTL_timestamps_Nlx(idx_good_TTL_Nlx))/1000; %if positive - it means that Nlx is ahead of neurologger

if isempty(idx_good_TTL_Nlg)
    disp('===========');
    disp(['ERROR in SYNC of Nlx TTLs for ',file_name_Nlx_TTL])
    disp('===========');
    error_sync=1;
    save(filename_error,'error_sync');
    return;
end;

%% Plotting the fit and the residuals (initial fit and final fit)

figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[21 2 25 15]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[-20.9 -0.6 0 0]);

% Plotting the initial fit
subplot(2,2,1)
hold on;
plot(TTL_timestamps_Nlg(idx_matching_Nlg_TTL)*1000,TTL_timestamps_Nlx(idx_matching_Nlx_TTL),'ob');
plot ([TTL_timestamps_Nlg(idx_matching_Nlg_TTL(1))*1000, TTL_timestamps_Nlg(idx_matching_Nlg_TTL(end))*1000], [polyval(initial_polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg(idx_matching_Nlg_TTL(1))*1000) , polyval(initial_polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg(idx_matching_Nlg_TTL(end))*1000)], '-k');
xlabel('Nlg timestamp (micro s)');
ylabel('Nlx timestamp (micro s)');
legend('Nlx vs. Nlg matching TTLs','Initial Linear Fit','Location','SouthEast');
title(sprintf('Initial sync of Nlx and Nlg clocks based on initial fit \n %d/%d  matching TTLs with maximal jitter < %d ms \n',length(idx_matching_Nlg_TTL),length(TTL_timestamps_Nlg), max_jitter_between_TTL_pairs));

%Plotting the residuals for the intial polyfit (calculated before removing the bad TTLs)
subplot(2,2,3)
hold on;
plot(1:length(matching_TTL_residuals_times),matching_TTL_residuals_times,'ob');
plot([1,length(matching_TTL_residuals_times)],[max_residual_TTL_fit_ms max_residual_TTL_fit_ms],'-r');
plot([1,length(matching_TTL_residuals_times)],[-max_residual_TTL_fit_ms -max_residual_TTL_fit_ms],'-r');
xlabel('# of matching TTLs');
ylabel('Time difference (ms)');
title(sprintf('Residuals of the intial fit \n (i.e. how much each actual measured TTL was different from the fit))'));
text(length(matching_TTL_residuals_times)/2,max_residual_TTL_fit_ms*0.8,'Cutoff for Good TTLs','Color',[1 0 0],'HorizontalAlignment','Center');

% Plotting the final fit
subplot(2,2,2)
hold on;
plot(TTL_timestamps_Nlg(idx_good_TTL_Nlg)*1000,TTL_timestamps_Nlx(idx_good_TTL_Nlx),'ob');
plot ([TTL_timestamps_Nlg(idx_good_TTL_Nlg(1))*1000, TTL_timestamps_Nlg(idx_good_TTL_Nlg(end))*1000], [polyval(polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg(idx_good_TTL_Nlg(1))*1000) , polyval(polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg(idx_good_TTL_Nlg(end))*1000)], '-k');
xlabel('Nlg timestamp (micro s)');
ylabel('Nlx timestamp (micro s)');
legend('Nlx vs. Nlg good TTLs','Final Linear Fit','Location','SouthEast');
title(sprintf('Final Sync of Nlx and Nlg clocks based on final (good) fit \n %d/%d good TTLs with residual < %d ms \n',length(idx_good_TTL_Nlg),length(TTL_timestamps_Nlg), max_residual_TTL_fit_ms));

%Plotting the residuals for the final polyfit (calculated before removing the bad TTLs)
subplot(2,2,4)
hold on;
plot(1:length(good_TTL_residuals_times),good_TTL_residuals_times,'ob');
xlabel('# of good TTLs');
ylabel('Time difference (ms)');
title(sprintf('Residuals of the final fit \n'));


% Saving the Figure
filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\3_Nlg_Nlx_clocks_sync_and_Residuals'];
%save the fig as .tiff
print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
% %save the fig as .fig
saveas(gca,filename_fig_to_save,'fig');


%% 6) For DBG: Computing the 'PSTH' for all detected TTLs and separately for the bad and the good TTL pulses, aligned to the detection point

PSTH_window_before_onset= round(0.1*CSC_Sampling_Rate_Hz); %0.1 sec before
PSTH_window_after_onset= round(0.2*CSC_Sampling_Rate_Hz); % %0.5 sec before
% FOR ALL TTLs
for ii_TTL=1:1:length(TTL_timestamps_Nlx)
    TTL_PSTH(:,ii_TTL) = Samples_reshaped( TTL_onset_idx(ii_TTL)-PSTH_window_before_onset : TTL_onset_idx(ii_TTL)+PSTH_window_after_onset );
end;
% FOR GOOD TTLs (based on the residuals)
counter=1;
for ii_TTL=idx_good_TTL_Nlx
    TTL_PSTH_good(:,counter) = Samples_reshaped( TTL_onset_idx(ii_TTL)-PSTH_window_before_onset : TTL_onset_idx(ii_TTL)+PSTH_window_after_onset );
    counter=counter+1;
end;
% FOR BAD TTLs (based on the residuals)
counter=1;
TTL_PSTH_bad=[];
for ii_TTL=idx_bad_TTL_Nlx
    TTL_PSTH_bad(:,counter) = Samples_reshaped( TTL_onset_idx(ii_TTL)-PSTH_window_before_onset : TTL_onset_idx(ii_TTL)+PSTH_window_after_onset );
    counter=counter+1;
end;

% Plotting the 'PSTH' of ALL TTL pulses, aligned to the detection point
figure;
subplot(3,3,2);
hold on;
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),TTL_PSTH,'-k'); %PSTH (all TTL events)
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),mean(TTL_PSTH,2),'-b'); % average TTL shape
plot(1000*[ PSTH_window_before_onset+1 PSTH_window_before_onset+1]./CSC_Sampling_Rate_Hz,[min(min(TTL_PSTH)) max(max(TTL_PSTH))],'-r');
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('All detected TTLs \n'));

% Plotting the 'PSTH' of the GOOD TTL pulses, aligned to the detection point
subplot(3,3,5);
hold on;
plot((1000*[1:length(TTL_PSTH_good)]./CSC_Sampling_Rate_Hz),TTL_PSTH_good,'-k'); %PSTH 
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),mean(TTL_PSTH,2),'-b'); % average TTL shape
plot(1000*[ PSTH_window_before_onset+1 PSTH_window_before_onset+1]./CSC_Sampling_Rate_Hz,[min(min(TTL_PSTH)) max(max(TTL_PSTH))],'-r');
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('Only Good TTLs (based on residuals) \n'));

% Plotting the 'PSTH' of the BAD TTL pulses, aligned to the detection point
subplot(3,3,8);
hold on;
if ~isempty(TTL_PSTH_bad)
    plot((1000*[1:length(TTL_PSTH_bad)]./CSC_Sampling_Rate_Hz),TTL_PSTH_bad,'-k'); 
end;
plot((1000*[1:length(TTL_PSTH)]./CSC_Sampling_Rate_Hz),mean(TTL_PSTH,2),'-b'); % average TTL shape
plot(1000*[ PSTH_window_before_onset+1 PSTH_window_before_onset+1]./CSC_Sampling_Rate_Hz,[min(min(TTL_PSTH)) max(max(TTL_PSTH))],'-r');
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('Only Bad TTLs (based on residuals) \n'));

% Saving the Figure
filename_fig_to_save=[folder_name '\Nlx_VT_and_TTL\4_TTL_PSTH_aligned_to_peak_detected'];
%save the fig as .tiff
print(filename_fig_to_save,['-f' num2str(gcf)],'-dtiff','-cmyk','-r300' );
% %save the fig as .fig
saveas(gca,filename_fig_to_save,'fig');

close all;

%% 7) Saving the fit

% Saving the fit and the TTLs
TTLs_and_Nlx2Nlg_clock_fits.TTL_timestamps_Nlg = TTL_timestamps_Nlg;
TTLs_and_Nlx2Nlg_clock_fits.TTL_timestamps_Nlx = TTL_timestamps_Nlx;
TTLs_and_Nlx2Nlg_clock_fits.idx_good_TTL_Nlg = idx_good_TTL_Nlg;
TTLs_and_Nlx2Nlg_clock_fits.idx_good_TTL_Nlx = idx_good_TTL_Nlx;
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlg2Nlx_microsec = polyfit_Nlg2Nlx_microsec;
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlx2Nlg_microsec = polyfit_Nlx2Nlg_microsec;

save(filename_data_to_save,'-struct','TTLs_and_Nlx2Nlg_clock_fits');



