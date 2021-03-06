function [error_sync]= NLG_PRE_sync_Nlx2Nlg_by_TTL_2_0(file_name_Nlx_TTL,p,folder_name)

%% Parameters
run_in_chunks_to_save_memory=0; % 1 - runs cross-corr in chunks to avoid memory overload or 0 - run in one piece (for stronger computers)
Nlg_event_delay_ms = 0; % Changed to '0' becasue TTL is now generated by transiever and there is no PC delay (used to be 35msec - This is the average delay introduced by the PC (between the emission of TTL pulse and the registration of the TTL-out event). 
                        % Positive value - means that the TTL pulse was emitted before the TTL-out event was saved in the log-file.
heaviside_template_window = 50; %in number of samples. Thi is the half window size of the step function (before the step, and the same length after the step).
cross_corr_detection_peak_shift = 50; %in units of samples (for heaviside template of length 100(50 -1 and 50 +1);
% cross_corr_threshold=6*10^5; %changed by Gily due to detection on non-TTLs %5*10^5; %used to identify the peaks of the crosscorrelation, corresponding to the rising phase of the TTL
cross_corr_threshold = 1*10^5; % changed by Maya due to non-detection on of TTLs %5*10^5; 
                             % used to identify the peaks of the crosscorrelation, corresponding to the rising phase of the TTL
max_jitter_between_TTL_pairs = 80;  %in ms  - this is the maximal acceptablle jitter between the Nlx and Nlg TTLs pairs, that will be used for synchronization
max_residual_TTL_fit_ms=15; % if the residual time between the data and the fit is bigger than this, we will recalculate the fit without this TTL time
filename_data_to_save=[folder_name '\Nlx_VT_and_TTL\TTLs_and_Nlx2Nlg_clock_fits.mat'];
filename_error=[folder_name '\Nlx_VT_and_TTL\error_sync.mat'];
TTL_tresh = 5;

%intializing.
error_sync=0;  %If there will be an error in sync (i.e. no good TTLs pairs) it will change to 1
TTL_timestamps_Nlx=[];
polyfit_Nlg2Nlx_microsec=[];
polyfit_Nlx2Nlg_microsec=[];

%% ======== Extract CSC min-max data  ============

% skip NLX TTL if already done
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


% sysntax: 
% [TTL_timestamps_Nlx, polyfit_Nlg2Nlx_microsec, polyfit_Nlx2Nlg_microsec, error_sync]= NLG_PRE_sync_Nlx2Nlg_by_TTL_2_0(file_name_Nlx_TTL,TTL_timestamps_Nlg,folder_name)
%
% Extracting TTL onset times fro Nlx, by cross-correlation-threshold
% detection of the continuously sampled TTL channel ('Audio2')
%
% David Omer Aug 2015, based on Arseny Finkestein's code
% 		bugs fixed:
%       	1. smoothing filter was normalized. 
%       	2. zero-padding signals before smoothing.
%       	3. the code accepts shorter Nlx TTL as well. 
%       	4. scaling before fitting to increase accuracy. 
%       	5. TTL detection using convolution. 
%           6. signal is z-scored to use STDs for thresholds 
%           
% 


%% Parameter
%dbstop if error
set(0,'DefaultFigureWindowStyle','normal');

Nlg_event_delay_ms=0;               % This is the average delay introduced by the PC (between the emission of TTL pulse and the registration of the TTL-out event). Positve value - means that the TTL pulse was emitted before the TTL-out event was saved in the log-file
max_jitter_between_TTL_pairs = 80;  %in ms  - this is the maximal acceptablle jitter between the Nlx and Nlg TTLs pairs, that will be used for synchronization
max_residual_TTL_fit_ms = 15;       % if the residual time between the data and the fit is bigger than this, we will recalculate the fit without this TTL time
folder_name = [p.path_dataout,p.path_day_dir,'\Preprocessing_SyncNlg&Nlx\'];

if ~exist(folder_name,'dir')
    mkdir(folder_name)
end

filename_data_to_save=[folder_name,'TTLs_and_Nlx2Nlg_clock_fits.mat'];
TTL_tresh = 5; 

%intializing.
error_sync=0;  %If there will be an error in sync (i.e. no good TTLs pairs) it will change to 1

%% 1) Extracting the TTLs from the Neurolynx CSC and from the p structure for the neurlogger

FieldSelection(1) = 1;
FieldSelection(2) = 1;
FieldSelection(3) = 1;
FieldSelection(4) = 1;
FieldSelection(5) = 1;
ExtractHeader = 1; %We want to extract the header

% we want to extract all 
ExtractMode = 1; % we want to extract the whole range
[Timestamps,~, SampleFrequency,~, Samples,~] = Nlx2MatCSC(file_name_Nlx_TTL, FieldSelection, ExtractHeader, ExtractMode) ;

CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
timev = 0:size(Samples,1)-1;
timev = repmat(timev(:),1,size(Samples,2));
timev = timev.*CSC_SamplePeriod_microsec;
Timestamps = repmat(Timestamps,size(Samples,1),1)+timev; 
Samples = Samples(:); 
Timestamps = Timestamps(:); 
Samples  = Samples- mean(Samples);  

%  detect ttl TTL events in the neuralynx data
% I switch from cross correlation to simple findpeak - detection of events with no offset. ; Didi 25.7.2015 
z = zscore(Samples(:));
heaviside_template_window = 50;
heaviside_template=[-ones(1,heaviside_template_window),ones(1,heaviside_template_window)]; %heaviside step_function imitating the rise time of the TTL pulse

c =conv(z,fliplr(heaviside_template),'same'); 
c = zscore(c(:)); 
[~,TTL_onset_idx] = findpeaks(c,'MinPeakHeight',TTL_tresh);
TTL_onset_idx =sort(TTL_onset_idx(:));
TTL_timestamps_Nlx = Timestamps(TTL_onset_idx);  % in microsec 
TTL_timestamps_Nlx = TTL_timestamps_Nlx(:)';

%%


% Extract Nlg TTL timestamps
TTL_timestamps_Nlg = (p.TTL_timestamps_Nlg'/1000)-Nlg_event_delay_ms; %in millisec

%% 2) Detecting missalignments (misses or false positive) between Nlx and Nlg TTLs
%Creating a vector from the first to the last TTL,  with bin size of 1 ms,

TTL_Nlx = zeros(1,ceil((TTL_timestamps_Nlx(end) - TTL_timestamps_Nlx(1))/1000)); 
TTL_Nlg = zeros(1,ceil(TTL_timestamps_Nlg(end)-TTL_timestamps_Nlg(1))); 

% Placing 1 at a timestamp (with bin of 1 ms) where a TTL was sent (Nlg) or detected (Nlx)
TTL_Nlx (round( (TTL_timestamps_Nlx - TTL_timestamps_Nlx(1))/1000) +1 )=1;
TTL_Nlg (round( (TTL_timestamps_Nlg - TTL_timestamps_Nlg(1))) +1 )=1;
 
% normelizing the kernel ; didi 25.8.2015
% very important the the kernal duration will be up to 1/10 of the average jitter
% used in the Neurolloger to produe the TTLs  
kernel = hamming(100)/sum(hamming(100));  

% paddign the vectors with zeros to avoid edge effect on the fitlering ; didi 25.8.2015
TTL_Nlx = [ zeros(1,length(kernel)*2),TTL_Nlx, zeros(1,length(kernel)*2)]; 
TTL_Nlg = [ zeros(1,length(kernel)*2),TTL_Nlg, zeros(1,length(kernel)*2)]; 

% smoothing with Hamming window(almost gaussian...;)) ; didi 25.8.2015
TTL_Nlg_smoothed=filtfilt(kernel,1,TTL_Nlg); 
TTL_Nlx_smoothed=filtfilt(kernel,1,TTL_Nlx);
TTL_Nlx_smoothed = TTL_Nlx_smoothed(length(kernel)*2 +1:end - (length(kernel)*2));  
TTL_Nlg_smoothed = TTL_Nlg_smoothed(length(kernel)*2 +1:end - (length(kernel)*2)); 
[c_Nlg_Nlx, lag] = xcorr(TTL_Nlg_smoothed,TTL_Nlx_smoothed);

% finding the delay in ms
[~,shift_Nlg_Nlx_ms ]=max(c_Nlg_Nlx);
shift_Nlg_Nlx_ms   = lag(shift_Nlg_Nlx_ms)+ TTL_timestamps_Nlg(1); % absolute shift!  (Didi) 
TTL_timestamps_Nlx_shift=(TTL_timestamps_Nlx- TTL_timestamps_Nlx(1)) + shift_Nlg_Nlx_ms*1000;

% remove zero padding from the vectors 
TTL_Nlx = TTL_Nlx(length(kernel)*2 +1:end - (length(kernel)*2));  
TTL_Nlg = TTL_Nlg(length(kernel)*2 +1:end - (length(kernel)*2)); 
axis_TTL_Nlg = linspace(0,length(TTL_Nlg)-1,length(TTL_Nlg)); 
axis_TTL_Nlx = linspace(0,length(TTL_Nlx)-1,length(TTL_Nlx));
axis_TTL_Nlx = axis_TTL_Nlx+ (shift_Nlg_Nlx_ms-TTL_timestamps_Nlg(1))+1;
intersect_idx = find(axis_TTL_Nlg>=axis_TTL_Nlx(1)-100 & axis_TTL_Nlg<=axis_TTL_Nlx(end)+100);

% Plotting crosscorrelation between the TTL timestamps in Nlx and Nlg

figure;
subplot(2,1,1);
hold on;
plot(lag,c_Nlg_Nlx,'-k'); %cross correlation
plot([ 0 0 ],[min(c_Nlg_Nlx) max(c_Nlg_Nlx)],'-.r');
%legend('Nlg Nlx crosscorr','0 lag');
title(['Computed peak shift =',num2str(shift_Nlg_Nlx_ms-TTL_timestamps_Nlg(1)), 'ms']);
xlabel('Lag [ms]');
subplot(2,1,2)
hold on
 
[c, lag2] = xcorr(TTL_Nlg_smoothed(intersect_idx),TTL_Nlx_smoothed);
plot(lag2,c,'-k'); %cross correlation
plot([ 0 0 ],[min(c_Nlg_Nlx) max(c_Nlg_Nlx)],'-.r');
xlabel('Lag [ms]');
title('cross-correlation, shift corrected'); 

% % Saving the Figure
filename_fig_to_save=[folder_name,'Nlg_Nlx_crosscorr'];
h = gcf; 
print(filename_fig_to_save,['-f' num2str(h.Number)],'-dtiff','-cmyk','-r300' );
saveas(gca,filename_fig_to_save,'fig');

%% 3) Aligning the TTLs (in case there were misses) and taking only the TTLs that were present both in Nlg and Nlx
% Updating the Nlx timestamps by the shift (in micro s), which could result
% from missalignment of the TTL trains i.e. - false positive/negative)

% detecting legitimate TTL pairs from Nlg and Nlx 
TTL_timestamps_Nlg = TTL_timestamps_Nlg(:)';
TTL_timestamps_Nlx_shift = TTL_timestamps_Nlx_shift(:)';
d  = abs(repmat(TTL_timestamps_Nlx_shift(:),1,length(TTL_timestamps_Nlg)) - repmat(TTL_timestamps_Nlg*1000,length(TTL_timestamps_Nlx_shift),1));  
[d,ii] = min(d');
ii2 = find(d/1000<=max_jitter_between_TTL_pairs);
d =d(ii2);
ii = ii(ii2); 
TTL_timestamps_Nlg_valid = TTL_timestamps_Nlg(ii);

%% 4) Computing the initial fit (for conversion of Nlx timestamps into Nlg timestamps) and residuals on matching TTLs
[initial_polyfit_Nlg2Nlx_microsec,~,muNlg2Nlx1] = polyfit(TTL_timestamps_Nlg_valid*1000, TTL_timestamps_Nlx , 1);
[initial_polyfit_Nlx2Nlg_microsec,~,muNlx2Nlg1] = polyfit(TTL_timestamps_Nlx , TTL_timestamps_Nlg_valid*1000,  1);

% Computing the residuals on matching TTLs ( i.e. how much each actual measured TTL was different from the fit)
matching_TTL_residuals_times =TTL_timestamps_Nlg_valid-  polyval(initial_polyfit_Nlx2Nlg_microsec,TTL_timestamps_Nlx,[],muNlx2Nlg1)/1000; %if positive - it means that Nlx is ahead of neurologger

% Identifying 'Good' and 'Bad' TTL pairs based on the residuals of the inital fit
%good
idx_good_TTL_Nlg= find(abs(matching_TTL_residuals_times)<max_residual_TTL_fit_ms);
idx_good_TTL_Nlx= find(abs(matching_TTL_residuals_times)<max_residual_TTL_fit_ms);
%bad
idx_bad_TTL_Nlg= find(abs(matching_TTL_residuals_times)>max_residual_TTL_fit_ms);
idx_bad_TTL_Nlx= find(abs(matching_TTL_residuals_times)>max_residual_TTL_fit_ms);

idx_valid_TTL_Nlg = ii(idx_good_TTL_Nlg); 
idx_valid_TTL_Nlx = idx_good_TTL_Nlx;



%% 5) Recalculating the final fit (for conversion of Nlx timestamps into Nlg timestamps) and residuals on Good TTL, i.e. - those TTL whose residuals didn't deviate significantly from the intial fit
[polyfit_Nlg2Nlx_microsec,~,muNlg2Nlx] = polyfit(TTL_timestamps_Nlg_valid(idx_good_TTL_Nlg)*1000, TTL_timestamps_Nlx(idx_good_TTL_Nlx), 1);
[polyfit_Nlx2Nlg_microsec,~,muNlx2Nlg] = polyfit(TTL_timestamps_Nlx(idx_good_TTL_Nlx),TTL_timestamps_Nlg_valid(idx_good_TTL_Nlg)*1000,  1);
clock_lag_Nlg2Nlx_micro_sec = polyval(polyfit_Nlg2Nlx_microsec,0,[],muNlg2Nlx); %if positive - it means that Nlx is ahead of neurologger

% Computing the residuals on final TTL - i.e. based on Good TTL
good_TTL_residuals_times =TTL_timestamps_Nlg_valid(idx_good_TTL_Nlg)-  polyval(polyfit_Nlx2Nlg_microsec,TTL_timestamps_Nlx(idx_good_TTL_Nlx),[],muNlx2Nlg)/1000; %if positive - it means that Nlx is ahead of neurologger

if isempty(idx_good_TTL_Nlg)
    disp('===========');
    disp(['ERROR in SYNC of Nlx TTLs for ',file_name_Nlx_TTL])
    disp('===========');
    error_sync=1;
    return;
end;

% Plotting the fit and the residuals (initial fit and final fit)

%%
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
plot(TTL_timestamps_Nlg_valid*1000,TTL_timestamps_Nlx,'ob');
plot ([TTL_timestamps_Nlg_valid(1)*1000, TTL_timestamps_Nlg_valid(end)*1000], [polyval(initial_polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg_valid(1)*1000,[],muNlg2Nlx1) , polyval(initial_polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg_valid(end)*1000,[],muNlg2Nlx1)], '-k');
xlabel('Nlg timestamp (micro s)');
ylabel('Nlx timestamp (micro s)');
legend('Nlx vs. Nlg matching TTLs','Initial Linear Fit','Location','SouthEast');
title(sprintf('Initial sync of Nlx and Nlg clocks based on initial fit \n %d/%d  matching TTLs with maximal jitter < %d ms \n',length(TTL_timestamps_Nlg_valid),length(TTL_timestamps_Nlg_valid), max_jitter_between_TTL_pairs));

% Plotting the residuals for the intial polyfit (calculated before removing the bad TTLs)
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
plot(TTL_timestamps_Nlg_valid*1000,TTL_timestamps_Nlx,'ob');
plot ([TTL_timestamps_Nlg_valid(1)*1000, TTL_timestamps_Nlg_valid(end)*1000], [polyval(polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg_valid(1)*1000,[],muNlg2Nlx) , polyval(polyfit_Nlg2Nlx_microsec,TTL_timestamps_Nlg_valid(end)*1000,[],muNlg2Nlx)], '-k');
xlabel('Nlg timestamp (micro s)');
ylabel('Nlx timestamp (micro s)');
legend('Nlx vs. Nlg good TTLs','Final Linear Fit','Location','SouthEast');
title(sprintf('Final Sync of Nlx and Nlg clocks based on final (good) fit \n %d/%d good TTLs with residual < %d ms \n',length(idx_good_TTL_Nlg),length(TTL_timestamps_Nlg), max_residual_TTL_fit_ms));

% Plotting the residuals for the final polyfit (calculated before removing the bad TTLs)
subplot(2,2,4)
hold on;
plot(1:length(good_TTL_residuals_times),good_TTL_residuals_times,'ob');
xlabel('# of good TTLs');
ylabel('Time difference (ms)');
title(sprintf('Residuals of the final fit \n'));


% % Saving the Figure
filename_fig_to_save=[folder_name 'Nlg_Nlx_clocks_sync_and_Residuals'];
h = gcf; 
print(filename_fig_to_save,['-f' num2str(h.Number)],'-dtiff','-cmyk','-r300' );
saveas(gca,filename_fig_to_save,'fig');


%% 6) For DBG: Computing the TTL  TTL onset-triggering for all detected TTLs and separately for the bad and the good TTL pulses


tt = -100*1000:300*1000;
xx = tt/1000; 
tt = round(tt/CSC_SamplePeriod_microsec); 

tt = repmat(tt(:),1,length(TTL_onset_idx)); 
TTLs = repmat(TTL_onset_idx(:)',size(tt,1),1);
TTLs = TTLs +tt; 
TTLs = Samples(TTLs); 


% Plotting the 'PSTH' of ALL TTL pulses, aligned to the detection point
figure;
subplot(2,2,1);
hold on;
plot(xx,TTLs,'-k'); %PSTH (all TTL events)
plot(xx,mean(TTLs,2),'-b'); % average TTL shape
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('All detected TTLs \n'));
set(gca,'xlim',[xx(1),xx(end)]); 


% Plotting the 'PSTH' of the GOOD TTL pulses, aligned to the detection point
subplot(2,2,3);
hold on;
plot(xx,TTLs(:,idx_good_TTL_Nlx),'-k'); %PSTH 
plot(xx,mean(TTLs(:,idx_good_TTL_Nlx),2),'-b'); % average TTL shape
xlabel('Time (ms)');
ylabel('CSC Amplitude');
title(sprintf('Only Good TTLs (based on residuals) \n'));
set(gca,'xlim',[xx(1),xx(end)]); 

% Plotting the 'PSTH' of the BAD TTL pulses, aligned to the detection point
subplot(2,2,4);
hold on;
if ~isempty(idx_bad_TTL_Nlx)
    plot(xx,TTLs(:,idx_bad_TTL_Nlx),'-k');
    
    plot(xx,mean(TTLs(:,idx_bad_TTL_Nlx),2),'-b'); % average TTL shape
    xlabel('Time (ms)');
    ylabel('CSC Amplitude');
    title(sprintf('Only Bad TTLs (based on residuals) \n'));
    set(gca,'xlim',[xx(1),xx(end)]); 
end;

% % Saving the Figure
filename_fig_to_save=[folder_name 'TTLs_aligned_to_peak_detected'];
h = gcf; 
print(filename_fig_to_save,['-f' num2str(h.Number)],'-dtiff','-cmyk','-r300' );
saveas(gca,filename_fig_to_save,'fig');


%% 7) Saving the fit

% Saving the fit and the TTLs
TTLs_and_Nlx2Nlg_clock_fits.TTL_timestamps_Nlg = TTL_timestamps_Nlg;
TTLs_and_Nlx2Nlg_clock_fits.TTL_timestamps_Nlx = TTL_timestamps_Nlx;
TTLs_and_Nlx2Nlg_clock_fits.idx_good_TTL_Nlg = idx_valid_TTL_Nlg;
TTLs_and_Nlx2Nlg_clock_fits.idx_good_TTL_Nlx = idx_valid_TTL_Nlx;
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlg2Nlx_microsec = polyfit_Nlg2Nlx_microsec;
TTLs_and_Nlx2Nlg_clock_fits.polyfit_Nlx2Nlg_microsec = polyfit_Nlx2Nlg_microsec;
TTLs_and_Nlx2Nlg_clock_fits.muNlx2Nlg = muNlx2Nlg; 
TTLs_and_Nlx2Nlg_clock_fits.muNlg2Nlx = muNlg2Nlx; 
 
save(filename_data_to_save,'-struct','TTLs_and_Nlx2Nlg_clock_fits');
set(0,'DefaultFigureWindowStyle','docked');


