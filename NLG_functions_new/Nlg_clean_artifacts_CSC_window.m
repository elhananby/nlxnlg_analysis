function Nlg_clean_artifacts_CSC_window(filename_in, filename_out, file_num, ii_file)

%---------------------------
% Maya G. (adapted from Michael Yartsev telemetry code and Arseny Finkelstein NLG code) 06/2015
% detect and save artifact on ONE CSC file
%--------------------------


%-----------------------------------------------------------------------
% PARAMETERS:
PLT = 0; % '1' if you we want to plot for DBG or '0' otherwise.
AWidth = 33; % Minimal interval (in samples) between two adjacent artifact-threshold-crossings:
% two threshold-crossings with smaller interval will be MERGED
% into ONE artifact. The reason to choose this number
% (equvalent to 2 ms) is that spikes close than 1 ms to
% artifact start or end are removed anyway later on.
min_artifact_length = 4; % Minimum length of artifact (assuming peak of spikes) -- spikes are ALWAYS SHORTER than this.
buff = 33; % buffer around artifacts ~ 1ms
Confidence_boundary_addition = 130; % Artifact threshold - as percentage of the maximal and minimal values in the crudly sorted sleep data
pos_threshold=1500; % positive threshold for artifact removal (in uv)
neg_threshold=-1500; % positive threshold for artifact removal (in uv)


% total artifact time in the session (in samples)
session_artifact_time = [0,0];


% total artifact time in this file (in samples), used in the original Telemetry code to find and select the quieter antenna
file_artifact_time = 0;

% load the current data file
current_filename	= [filename_in,num2str(ii_file)];
filename_out_save	= [filename_out,num2str(ii_file)];
if isempty(dir(sprintf('%s.mat',current_filename)))
    warning('file %d is missing',ii_file)
    return
end

eval(['load ', current_filename]);		  % Load the current CSC data
current_samples = CSC.data.Samples_filtered_saved;
current_timestamps = CSC.data.Timestamps_filtered_samples_saved;
clear CSC

Artifacts_IX = {};
Artifacts_start_end_IX = {};
Artifacts_start_end_IX_no_removal = {};
Artifacts_start_end_timestamps = {};
Artifacts_start_end_timestamps_no_removal = {};

for ii_channel = 1:4
    
    disp(['Detecting artifacts on channel #',num2str(ii_channel),...
        ' of file #', num2str(ii_file),							...
        ' out of ',num2str(file_num),                           ...
        ' files...']);
    
    if (isempty(current_samples))
        continue
    end
    
    Samples = current_samples{ii_channel};
    
    pos_artifact_points = find( Samples > pos_threshold);			% get positive artifacts
    neg_artifact_points = find( Samples < neg_threshold);			% get negative artifacts
    artifact_points = sort(unique([pos_artifact_points neg_artifact_points]));
    
    file_artifact_time = file_artifact_time + length(artifact_points);
    
    if ~isempty(artifact_points)	% Check  if we had artifacts at all
        
        % plot debug info
        if PLT
            h = figure;
            set(gcf,'DefaultAxesFontSize',8);
            set(gcf,'DefaultAxesFontName','helvetica');
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 22 21]);
            set(gcf,'PaperOrientation','portrait');
            set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[4 5 0 0]);
            a1 = subplot(2,1,1);
            plot(Samples,'k')
            ylabel('signal of this flight')
            title('Before Artifact Removal')
            hold all
            
            Samples_no_wing = Samples;	% for later debug plot
        end
        
        artifact_counter = 0;
        merged_artifact_points = [];
        
        % go over the artifact points, merging close points
        for curr_artifact_point_idx = 1:length(artifact_points);
            curr_artifact_point = artifact_points(curr_artifact_point_idx);
            merged_artifact_points = [merged_artifact_points curr_artifact_point];
            
            % check if no additional points need to be merged into this artifact
            if  (curr_artifact_point == artifact_points(end) || curr_artifact_point + AWidth <  artifact_points(curr_artifact_point_idx + 1))
                % check if the artifact is too short, if so ignore it
                if length(merged_artifact_points) < min_artifact_length
                    merged_artifact_points = [];									% reinitialize merged artifact
                else
                    artifact_counter = artifact_counter + 1;
                    first_artifact_point = merged_artifact_points(1);
                    last_artifact_point = merged_artifact_points(end);
                    
                    Zero_cross_start_without_removal = max(find(Samples(1:first_artifact_point)*sign(Samples(first_artifact_point)) < 0));
                    Zero_cross_end_without_removal	 = last_artifact_point + min(find(Samples(last_artifact_point:end)*sign(Samples(last_artifact_point)) < 0));
                    Zero_cross_start				 = Zero_cross_start_without_removal - buff;
                    Zero_cross_end					 = Zero_cross_end_without_removal   + buff;
                    
                    if Zero_cross_end > length(current_timestamps)
                        Zero_cross_end = length(current_timestamps);
                    end
                    
                    if Zero_cross_start < 1
                        Zero_cross_start = 1;
                    end
                    
                    Artifacts_IX{ii_channel,artifact_counter} = [Zero_cross_start:Zero_cross_end];
                    Artifacts_start_end_IX{ii_channel,artifact_counter} = [Zero_cross_start,Zero_cross_end];
                    Artifacts_start_end_timestamps{ii_channel,artifact_counter} = current_timestamps(1,[Zero_cross_start,Zero_cross_end]);
                    Artifacts_IX__no_removal{ii_channel,artifact_counter} = [Zero_cross_start_without_removal:Zero_cross_end_without_removal];
                    Artifacts_start_end_IX_no_removal{ii_channel,artifact_counter} = [Zero_cross_start_without_removal,Zero_cross_end_without_removal];
                    if Zero_cross_end_without_removal>length(current_timestamps) % Arseny (I added it to avoid crush)
                        Artifacts_start_end_timestamps_no_removal{ii_channel,artifact_counter} = current_timestamps(1,[Zero_cross_start_without_removal,length(current_timestamps)]);
                    else
                        Artifacts_start_end_timestamps_no_removal{ii_channel,artifact_counter} = current_timestamps(1,[Zero_cross_start_without_removal,Zero_cross_end_without_removal]);
                    end
                    merged_artifact_points = [];											% reinitialize merged artifact
                    
                    % plot debug info
                    if PLT
                        subplot(2,1,1)
                        hold on
                        plot(Artifacts_IX{ii_channel,artifact_counter},Samples(Artifacts_IX{ii_channel,artifact_counter}),'.:r')
                        axis([-inf inf min(Samples)-5 max(Samples)+5])
                        Samples_no_wing(Artifacts_IX{ii_channel,artifact_counter}) = 0;		% replace with zero's for later debug visualization purposes
                    end
                end
            end
        end
        
        % plot debug info
        if PLT
            a2 = subplot(2,1,2);
            plot(Samples_no_wing,'k')
            axis([-inf inf min(Samples)-5 max(Samples)+5])
            title('After Artifact Removal')
            linkaxes([a1 a2], 'x');
        end
    else												% If there were no artifacts
        Artifacts_start_end_IX{ii_channel,1} = [];
    end
    
end

Artifacts.params.pos_threshold			 = pos_threshold;
Artifacts.params.neg_threshold			 = neg_threshold;
Artifacts.params.AWidth					 = AWidth;
Artifacts.params.min_artifact_length	 = min_artifact_length;
Artifacts.params.buff					 = buff;

Artifacts.params.associated_file_name_in = current_filename;

Artifacts.data.Artifacts_start_end_IX					 = Artifacts_start_end_IX;
Artifacts.data.Artifacts_start_end_timestamps			 = Artifacts_start_end_timestamps;
Artifacts.data.Artifacts_start_end_IX_no_removal		 = Artifacts_start_end_IX_no_removal;
Artifacts.data.Artifacts_start_end_timestamps_no_removal = Artifacts_start_end_timestamps_no_removal;
Artifacts.data.Total_artifact_time = file_artifact_time;

eval( ['save ', filename_out_save, ' Artifacts'] );
fprintf('Total artifact time in minute %d : %d\n', ii_file, file_artifact_time)

end % func

