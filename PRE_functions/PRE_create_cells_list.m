function cell_files_to_load = PRE_create_cells_list(excel_file, P)
%PRE_CREATE_CELLS_LIST Automatically creates list of cells from sorted
%spike file
%   INPUT:  excel_file - excel file to open
%           P - structure with ntt files to open
%   OUTPUT: updated excel file
%           exp_files_to_load - list of experiment files to load, which
%           contain all the cell data needed for analysis

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
    
    % initiate temporary p structure
    p = P(ii_exp);
    
    % find and load ntt files of exp
    file_search = sprintf('%s\\%s\\*.ntt', p.path_dataout, p.datadir_out);
    files_to_load = subdir(file_search);
    
    % load video data file of exp
    vt_file_search = fullfile(p.path_dataout, p.datadir_out, 'VT*.mat');
    vt_files_to_load = subdir(vt_file_search);
    load(vt_files_to_load.name, 'vt');
    
    p.vtfile = vt_files_to_load.name;
    
    %% loop over sessions
    for ii_nses = 1:length(p.S)
        sessionStartTime = p.S(ii_nses).start_time;
        sessionEndTime = p.S(ii_nses).end_time;
        
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
            
            [Timestamps, CellNumbers, Samples] = Nlx2MatSpike( Filename, [1 0 1 0 1], 0, ExtractionMode, ExtractionModeVector); % load current spike file
            
            % if file not spike sorted
            if all(CellNumbers == 0)
                fprintf('No cells in file %s.\nPerform spike sorting and then press any key to continue\n', Filename);
                pause;
            end
            
            % if the cluster size = all spikes, this is a marker for a
            % spikeless file and we can continue
            if length(find(CellNumbers == 1)) == length(CellNumbers)
                continue;
            end
            
            cn = unique(CellNumbers(CellNumbers ~= 0)); % find all non-zero cell numbers
            
            %% loop over cells in CellNumbers
            for ii_cell = cn
                
                % create new cells line
                c = struct;
                
                c.cell_number = curr_cell;
                c.animal = p.animal;
                c.animal_name = p.animal_name;
                c.day = p.day;
                c.experiment = p.experiment;
                c.session = ii_nses;
                c.TT = ii_tetrode;
                c.cell_id = ii_cell;
                
                excel_line = table2cell(struct2table(c)); % create line to append to excel file
                
                %% interpolate spikes
                c = interp_spikes(c, Timestamps, CellNumbers, ii_cell, vt, p);
                
                % get only current session video info
                cVt = session_video(sessionStartTime, sessionEndTime, vt);
                
                % claculate spike train
                c = spike_train(c, Timestamps, CellNumbers, ii_cell, vt);
                
                % get spike shape
                c.spikeShape = Samples(:, :, CellNumbers == ii_cell);
                
                %% save cell file and update excel sheet
                
                % create index to check if cell already exists in db
                current_cell_string = sprintf('%i_%s_%i_%i_%i_%i_%i',...
                    p.animal,...
                    p.animal_name,...
                    p.day,...
                    p.experiment,...
                    ii_nses,...
                    ii_tetrode,...
                    ii_cell);
                
                % create filename structure
                outfile_name = sprintf('%i_%i-%s_%s_%d_Exp%i_Session%i_TT%i_Cell%i.mat',...
                    c.cell_number, p.animal, p.animal_name, p.nlgnlx, p.day, p.experiment, ii_nses, ii_tetrode, ii_cell);
                
                outfile_FULL = fullfile(p.path_dataout,...
                    'Cells',...
                    sprintf('%i_%s', p.animal, p.animal_name),...
                    outfile_name);
                
                % check if the folder exists
                if ~exist(fileparts(outfile_FULL), 'dir')
                    mkdir(fileparts(outfile_FULL));
                end
                
                % no excel AND yes file
                if exist(outfile_FULL, 'file') && ~any(strcmp(current_cell_string, existing_cells))
                    
                    % update excel but don't save file
                    xlswrite(excel_file, excel_line, 'Cells', sprintf('A%i', curr_cell+1));
                    
                    % no file AND yes excel
                elseif ~exist(outfile_FULL, 'file') && any(strcmp(current_cell_string, existing_cells))
                    
                    % save file but don't update excel
                    save(outfile_FULL,  'c', 'cVt' , 'p');
                    
                    % no file AND no excel
                elseif ~exist(outfile_FULL, 'file') && ~any(strcmp(current_cell_string, existing_cells))
                    
                    % save file and update excel
                    xlswrite(excel_file, excel_line, 'Cells', sprintf('A%i', curr_cell+1));
                    save(outfile_FULL,  'c', 'cVt' , 'p');
                    
                end
                
                curr_cell = curr_cell + 1; % for excel file
                
                cell_files_to_load{end+1} = outfile_FULL;
            end % cell
            
        end % tetrode
        
    end % session
    
end % experiment

end

function c = interp_spikes(c, ts, cn, cell, vt, p)
%% INTERPSPIKES interpolates spike values
%   c - cell structure to update
%   ts - spikes timestamps
%   cn - CellNumbers array
%   cell - current processed cell
%   vt - video data

c.timestamps = ts(cn == cell)';

if any(p.throw_away_times)
    c.timestamps = c.timestamps(c.timestamps < p.throw_away_times(1) | c.timestamps > p.throw_away_times(2));
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
end

function c = spike_train(c, ts, cn, cell, vt)
% Calculates the spiketrain
dt = mean(diff(vt.timestamps));
timebins = [vt.timestamps; (vt.timestamps(end) + dt)];
c.spikeTrain = histcounts(ts(cn == cell), timebins)';

% Smooths firing rate
filter = gaussmf(-4:4, [2 0]); filter = filter / sum(filter);
c.firingRate = c.spikeTrain / dt;
c.smoothFiringRate = conv(c.firingRate, filter, 'same');
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