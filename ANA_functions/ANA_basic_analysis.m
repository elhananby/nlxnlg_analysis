function ANA_basic_analysis(cell_list)
dbstop if error;
close all; clc;

global m n count dt nPosBins nHdBins nSpeedBins lag boxSize;
m = 3; n = 4; % subplots
nPosBins = 30; nHdBins = 60; nSpeedBins = 10; lag = 500;

%% enumerate over cells and load one cell at a time

for ii_cell = cell_list
    
    load(ii_cell{1});
    fprintf('Animal %i-%s\tDay %i\tExperiment %i\tSession %i\tTetrode %i\tCell %i\n',...
        p.animal, p.animal_name, p.day, p.experiment, c.session, c.TT, c.cell_id);
    
    % correct struct names for analysis
    c = c;
    vt = cVt;
    
    p = p;
    s = p.S(c.session);
    
    boxSize = p.arena_width_east_to_west;
    dt = mean(diff(vt.timestamps));
    
    fig = figure('Units', 'Normalized',...
        'Position', [0 0 1 1],...
        'Visible', 'off');
    
    %% cleaning thresholds
    vtKeepIdx = index_to_keep(vt, p, s);
    cKeepIdx = index_to_keep(c, p, s);
    
    %% calculate basic stuff
    % rate map
    [posOccupancy, posSpikes, posRates, posRatesSmooth] = calculate_rate_map(c, vt, cKeepIdx, vtKeepIdx);
    
    % hd map
    [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx);
    
    % speed map
    [speedOccupancy, speedSpikes, speedRates] = calculate_speed_map(c, vt, cKeepIdx, vtKeepIdx);
    
    % isi
    spikeISI = diff(c.timestamps.*1e-6); % in seconds
    
    % xcorr
    [rSpikeTrain, lagsSpikeTrain] = xcorr(c.spikeTrain, 500);
    
    % border score
    borderScore = calculate_border_score(posRates);
    
    %% plot all the basic stuff
    count = 1;
    
    PLOT_behavior;
    PLOT_rate_map;
    PLOT_spike_shape;
    PLOT_ISI;
    PLOT_hd_time;
    PLOT_hd_polar;
    PLOT_hd_speed;
    PLOT_hd_histogram;
    PLOT_hd_shuffle;
    PLOT_border_shuffle;
    PLOT_speed_map;
    
    %% title
    suptitle(sprintf('%s \t Animal %i \t Day %i \n Experiment %i \t Session %i - %s \n TT %i \t Cell %i',...
        p.nlgnlx,...
        p.animal,...
        p.day,...
        p.experiment,...
        c.session,...
        p.S(c.session).session,...
        c.TT,...
        c.cell_id));
    
    % save
    [filepath_fig, ~, ~] = fileparts(ii_cell{1});

    filename_fig = sprintf('%i_%i-%s_%s_Day%d_Exp%i_Session%i_TT%i_Cell%i',...
        c.cell_number, p.animal, p.animal_name, p.nlgnlx, p.day, p.experiment, c.session, c.TT, c.cell_id);
    
    savefig(fig, fullfile(filepath_fig, filename_fig), 'compact');
    saveas(fig, fullfile(filepath_fig, [filename_fig '.png']), 'png');
    close all;
    fprintf('\t\tSaved %s\n', fullfile(filepath_fig, filename_fig));
end % cells

end % functions
