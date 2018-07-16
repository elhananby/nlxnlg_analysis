function [  polyfit_Nlg2Nlx_microsec polyfit_Nlx2Nlg_microsec ]= NLG_PRE_sync_Nlx2Nlg_by_LED (folder_name,p)
% =======================================================================

%---------------------------
% Arseny Finkelstein 18/03/2014

% Synchronization of Neuralynx VT with Neurologger By LEDs on/off times (in case we don't have TTLs)

%--------------------------

filename_data_to_save=[folder_name '\Nlx_VT_and_TTL\TTLs_and_Nlx2Nlg_clock_fits.mat'];

% skip CSC Artifact cleaning if we had done this already
if exist ((filename_data_to_save))==2
    disp('======================');
    disp (['Skipping Synchronizing by LEDs. Already extracted in ' filename_data_to_save]);
    disp('======================');
    load(filename_data_to_save);
    return;
else
    % =======================================================================
    disp('===========');
    disp(['Synchronizing by LEDs for ',file_name_Nlx_TTL])
    disp('===========');
end;

%% 1)  Get Nlg and Nlx sync times 

%Getting Nlg sync times
disp('Showing NLG events for LED sync')
p.Nlg_EventStrings
Nlg_LED_turn_on_events=[5;17];
Nlg_LED_turn_off_events=[6;18];
Nlg_LED_turn_on_times=cell2mat(p.Nlg_EventTimestamps(Nlg_LED_turn_on_events))'*1000; %in microsec
Nlg_LED_turn_off_times=cell2mat(p.Nlg_EventTimestamps(Nlg_LED_turn_off_events))'*1000;%in microsec

%Getting Nlx sync times
disp('Showing NLG events for LED sync')
p.Nlx_EventStrings
Nlx_LED_turn_on_events=[2;13];
Nlx_LED_turn_off_events=[3;14];
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
    Nlx_LED_off_Timestamps_v1(ii_sync)=Timestamps_v1(find(ExtractedX,1,'last')+1);
    
    VT2_file = fullfile(VT_directory,'VT2.nvt');
    [Timestamps_v2, ExtractedX, ExtractedY, Targets_VT2, NlxHeader] = ...
        Nlx2MatVT(VT2_file, Extract_Fields, Extract_Header, Extraction_Mode,ExtractionModeArray ) ; % Extract data - camera 2
    
    Nlx_LED_on_Timestamps_v2(ii_sync)=Timestamps_v2(find(ExtractedX,1,'first'));
    Nlx_LED_off_Timestamps_v2(ii_sync)=Timestamps_v2(find(ExtractedX,1,'last')+1);
    
end;

Nlx_LED_on_Timestamps_min_v1_v2=min([Nlx_LED_on_Timestamps_v1;Nlx_LED_on_Timestamps_v2])
Nlx_LED_off_Timestamps_min_v1_v2=min([Nlx_LED_off_Timestamps_v1;Nlx_LED_off_Timestamps_v2])

LED_timestamps_Nlg = [Nlg_LED_turn_on_times, Nlg_LED_turn_off_times];
LED_timestamps_Nlx = [Nlx_LED_on_Timestamps_min_v1_v2, Nlx_LED_off_Timestamps_min_v1_v2];

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


