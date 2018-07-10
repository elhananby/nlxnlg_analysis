function cell_files_to_load = PRE_create_cells_list(excel_file, P)
%PRE_CREATE_CELLS_LIST Automatically creates list of cells from sorted
%spike file
%   INPUT:  excel_file - excel file to open
%           P - structure with ntt files to open
%   OUTPUT: updated excel file
%           exp_files_to_load - list of experiment files to load, which
%           contain all the cell data needed for analysis
dbstop if error

T = readtable(excel_file,...
    'Sheet', 2,...
    'ReadVariableNames', 1);

if isempty(T) % meaning, this is a completely empty worksheet
    existing_cells = {};
    curr_cell = 1;
else
    
    TC = table2cell(T);
    
    for ii = 1:height(T)
        existing_cells{ii} = sprintf('%i_%s_%i_%i_%i_%i_%i',...
            TC{ii, 2}, TC{ii, 3}, TC{ii, 4}, TC{ii, 5}, TC{ii, 6}, TC{ii, 7}, TC{ii, 8});
    end
    
    curr_cell = height(T) + 1;
end

cell_files_to_load = {};

%% loop over experiments
for ii_exp = 1:length(P)
    fprintf('Exp %i %.2f%%\n', ii_exp, ii_exp*100/length(P));
    % initiate temporary p structure
    p = P(ii_exp);
    
    % find and load ntt files of exp
    file_search = sprintf('%s\\%s\\*.ntt', p.path_dataout, p.datadir_out);
    files_to_load = subdir(file_search);
    
    % load video data file of exp
    vt_file_search = fullfile(p.path_dataout, p.datadir_out, 'VT*.mat');
    vt_files_to_load = subdir(vt_file_search);
    load(vt_files_to_load.name, 'vt');
    orgVt = vt;
    
    p.vtfile = vt_files_to_load.name;
    
    %% loop over sessions
    for ii_nses = 1:length(p.S)
        sessionStartTime = p.S(ii_nses).start_time;
        sessionEndTime = p.S(ii_nses).end_time;
        s = p.S(ii_nses);
        
        %% loop over tetrode files
        for ii_tetrode = 1:length(files_to_load)
            
            % load ntt
            Filename = files_to_load(ii_tetrode).name;
            if length(p.S) == 1 % only one sessions
                ExtractionMode = 1;
                ExtractionModeVector = [];
            else
                ExtractionMode = 4;
                ExtractionModeVector = [sessionStartTime sessionEndTime]; % load only the session time we want
            end
            
            [Timestamps, CellNumbers, FD, Samples, Header] = Nlx2MatSpike( Filename, [1 0 1 1 1], 1, ExtractionMode, ExtractionModeVector); % load current spike file
            
    
            
            % if file not spike sorted
            if all(CellNumbers == 0)
                continue;
%                 line = sprintf('No cells in %s.\nIs this intentional? Y/N [Y]. ', Filename);
%                 str = input(line, 's');                
%                 if isempty(str), str = 'y'; end
%                 fprintf(repmat('\b', 1, length(line)));
%                 if strcmpi(str, 'Y'), continue;
%                 elseif strcmpi(str, 'N'), fprintf('Perform Spike sorting and then press any key to continue'); pause; end
            end
            
            % if the cluster size = all spikes, this is a marker for a spikeless file and we can continue
            if length(find(CellNumbers == 1)) == length(CellNumbers)
                continue;
            end
            
            cn = unique(CellNumbers(CellNumbers ~= 0)); % find all non-zero cell numbers
            
            %% loop over cells in CellNumbers
            for ii_cell = cn
                
                %% check if cell exists in excel
                % ID, Name, Day, Exp, Session, TT, Cell
                % create index to check if cell already exists in db
                current_cell_string = sprintf('%i_%s_%i_%i_%i_%i_%i',...
                    p.animal,...        % ID
                    p.animal_name,...   % name
                    p.day,...           % day
                    p.experiment,...    % exp
                    ii_nses,...         % session
                    ii_tetrode,...      % TT
                    ii_cell);           % cell
                
                % check if line already exists in excel cell list
                curr_existing_cell = find(strcmp(current_cell_string, existing_cells));
                
                %% create new cells line
                metaData = struct;
                
                % if the line already exists in cell list, we want to take
                % the cell number from there
                if ~isempty(curr_existing_cell)
                    metaData.cell_number = curr_existing_cell;
                else % otherwise add new cell using curr_cell
                    metaData.cell_number = curr_cell;
                end
                
                metaData.animal = p.animal;
                metaData.animal_name = p.animal_name;
                metaData.day = p.day;
                metaData.experiment = p.experiment;
                metaData.session = ii_nses;
                metaData.TT = ii_tetrode;
                metaData.cell_id = ii_cell;  
                
                excel_line = table2cell(struct2table(metaData)); % create line to append to excel file
                
                %% cluster quality
                ClusterSpikes = find(CellNumbers == ii_cell);
                metaData.IsoDist = IsolationDistance(FD', ClusterSpikes);
                metaData.Lratio = L_Ratio(FD', ClusterSpikes);
                
                %% interpolate spikes
                c = interp_spikes(Timestamps, CellNumbers, ii_cell, orgVt, p);
                c = index_to_keep(c, p, s);
                
                % get only current session video info
                vt = session_video(sessionStartTime, sessionEndTime, orgVt);
                vt = index_to_keep(vt, p, s);
                
                % claculate spike train
                % c = spike_train(c, Timestamps, CellNumbers, ii_cell, vt);
                
                % get spike shape (if not nlg session)
                if ~strcmpi(p.nlgnlx, 'nlg')
                    
                    ADBits2Volts = regexpi(Header(contains(Header, 'ADBitVolts')),...
                        '(\d.\d+)','match'); % get adbits2volts from header
                    ADBits2Volts = cellfun(@str2double, ADBits2Volts{:}) * 1e6; % convert to double microvolt
                else % then the values are already in microvolt
                    
                    ADBits2Volts = [1 1 1 1];
                end
                
                spikeShape = [];
                for iiChan = 1:size(Samples, 2)
                    spikeShape(:, iiChan, :) = Samples(:, iiChan, CellNumbers == ii_cell) .* ADBits2Volts(iiChan); % convert to microvolts
                end
                
                %% save cell file and update excel sheet
                % create filename structure
                outfile_name = sprintf('%i_%i-%s_%s_Day%d_Exp%i_Session%i_TT%i_Cell%i.mat',...
                    metaData.cell_number, p.animal, p.animal_name, p.nlgnlx, p.day, p.experiment, ii_nses, ii_tetrode, ii_cell);
    
                outfile_FULL = fullfile(p.path_dataout,...
                    'Cells',...
                    sprintf('%i_%s', p.animal, p.animal_name),...
                    outfile_name);
                
                % check if the folder exists
                if ~exist(fileparts(outfile_FULL), 'dir')
                    mkdir(fileparts(outfile_FULL));
                end
                
                % yes file AND no excel
                if exist(outfile_FULL, 'file') && ~any(strcmp(current_cell_string, existing_cells))
                    
                    % update excel but don't save file
                    xlswrite(excel_file, excel_line, 'Cells', sprintf('A%i', curr_cell+1));
                    curr_cell = curr_cell + 1; % for excel file
                    
                % no file AND yes excel
                elseif ~exist(outfile_FULL, 'file') && any(strcmp(current_cell_string, existing_cells))
                    
                    % save file but don't update excel
                    save(outfile_FULL, 'metaData', 'c', 'vt', 'p', 'spikeShape');
                    
                % no file AND no excel
                elseif ~exist(outfile_FULL, 'file') && ~any(strcmp(current_cell_string, existing_cells))
                    
                    % save file and update excel
                    xlswrite(excel_file, excel_line, 'Cells', sprintf('A%i', curr_cell+1));
                    save(outfile_FULL, 'metaData', 'c', 'vt', 'p', 'spikeShape');
                    curr_cell = curr_cell + 1; % for excel file
                    
                end
                
                cell_files_to_load{end+1} = outfile_FULL;

            end % cell
            
        end % tetrode
        
    end % session
    
end % experiment

end

function c = interp_spikes(ts, cn, cell, vt, p)
dbstop if error
%% INTERPSPIKES interpolates spike values
%   c - cell structure to update
%   ts - spikes timestamps
%   cn - CellNumbers array
%   cell - current processed cell
%   vt - video data

if strcmpi(p.nlgnlx, 'nlg')
    c.timestamps = ts(cn == cell)' + polyval(p.nlg.align_timestamps_nlg2nlx.p,...
        ts(cn == cell)',...
        p.nlg.align_timestamps_nlg2nlx.S,...
        p.nlg.align_timestamps_nlg2nlx.mu);
else, c.timestamps = ts(cn == cell)';
end

c.posx       = interp1(vt.timestamps, vt.posx, c.timestamps);
c.posx2      = interp1(vt.timestamps, vt.posx2, c.timestamps);
c.posy       = interp1(vt.timestamps, vt.posy, c.timestamps);
c.posy2      = interp1(vt.timestamps, vt.posy2, c.timestamps);
c.posx_c     = interp1(vt.timestamps, vt.posx_c, c.timestamps);
c.posy_c     = interp1(vt.timestamps, vt.posy_c, c.timestamps);
c.poshd      = interp1(vt.timestamps, vt.poshd, c.timestamps);
c.vx         = interp1(vt.timestamps, vt.vx, c.timestamps);
c.vy         = interp1(vt.timestamps, vt.vy, c.timestamps);
c.speed      = interp1(vt.timestamps, vt.speed, c.timestamps);
c.angVel     = interp1(vt.timestamps, vt.angVel, c.timestamps);
end

function c = spike_train(c, ts, cn, cell, vt)
% Calculates the spiketrain
dt = mean(diff(vt.timestamps));
timebins = [vt.timestamps; (vt.timestamps(end) + dt)];
metaData.spikeTrain = histcounts(ts(cn == cell), timebins)';

% Smooths firing rate
filter = gaussmf(-4:4, [2 0]); filter = filter / sum(filter);
metaData.firingRate = metaData.spikeTrain / dt;
metaData.smoothFiringRate = conv(metaData.firingRate, filter, 'same');
end

function vt = session_video(startTime, endTime, vt)
% get only session part of video

sIdx = vt.timestamps >= startTime & vt.timestamps <= endTime;
vt.timestamps = vt.timestamps(sIdx);
vt.posx = vt.posx(sIdx);
vt.posy = vt.posy(sIdx);
vt.posx2 = vt.posx2(sIdx);
vt.posy2 = vt.posy2(sIdx);
vt.posx_c = vt.posx_c(sIdx);
vt.posy_c = vt.posy_c(sIdx);
vt.poshd = vt.poshd(sIdx);
vt.vx = vt.vx(sIdx);
vt.vy = vt.vy(sIdx);
vt.speed = vt.speed(sIdx);
end

function IsoDist = IsolationDistance(FD, ClusterSpikes)

% IsoDist = IsolationDistance(FD, ClusterSpikes)
%
% Isolation Distance
% Measure of cluster quality
%
% Inputs:   FD:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into FD which lists spikes from the cell whose quality is to be evaluated.
%
% Created by Ken Harris
% 
% Code by ADR 2012/12, from earlier versions

[nSpikes, nCh] = size(FD);

nClusterSpikes = length(ClusterSpikes);

if nClusterSpikes > nSpikes/2
    IsoDist = nan; % not enough out-of-cluster-spikes - IsoD undefined
    return
end

InClu = ClusterSpikes;
OutClu = setdiff(1:nSpikes, ClusterSpikes);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
m = mahal(FD, FD(ClusterSpikes,:));

mNoise = m(OutClu); % mahal dist of all other spikes

% calculate point where mD of other spikes = n of this cell
sorted = sort(mNoise);
IsoDist = sorted(nClusterSpikes);
end

function [output, m] = L_Ratio(FD, ClusterSpikes, NoiseSpikes)

% output = L_Ratio(FD, ClusterSpikes)
%
% L-ratio
% Measure of cluster quality
%
% Inputs:   FD:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into FD which lists spikes from the cell whose quality is to be evaluated.
%
% Output: a structure containing three components
%           Lratio, L, df
% ADR 2017-10-02 Added code to allow processing of NoiseSpikes

% find # of spikes in this cluster
[nSpikes, nD] = size(FD);

nClusterSpikes = length(ClusterSpikes);

% mark spikes which are not cluster members
if nargin < 3
    NoiseSpikes = setdiff(1:nSpikes, ClusterSpikes);
end

m = mahal(FD, FD(ClusterSpikes,:));
df = size(FD,2);

L = sum(1-chi2cdf(m(NoiseSpikes),df));
Lratio = L/nClusterSpikes;

output.L = L;
output.Lratio = Lratio;
output.df = nD;
end