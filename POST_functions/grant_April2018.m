dbstop if error
clearvars
close all
global parms
parms.sigma = 3; % gaussiam smoothing factor
parms.time_per_bin= 0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
this is now manually calculated in line 55 for each S variable.

Minimum radius used in the auto-correlogram when finding the best
parms.bin_size = 3; % size of spatial bin (for create the rate map)
parms.num_of_direction_bins = 60; % for head-direction calculation
parms.max_lag = 500; % max lag (in msec) for temporal autocorrelation
export_folder = 'D:\Dropbox (Technion Dropbox)\Lab_folder\Grant_2018_data\HD';

% Plotting for grant
s = tdfread('D:\hdcells.txt');
cells = s.cell_num(logical(s.use));
cells = [177; 413];

for ii = 1:length(cells)
    
    % find data file
    file_search = sprintf('D:\\experiment_data\\cells\\%s_*.mat', num2str(cells(ii)));
    file_to_load = subdir(file_search);
    
    % load cell file
    load(file_to_load.name);
    
    posx = C.S.pos.x;
    posy = C.S.pos.y;
    
    spkx = C.S.spk.x;
    spky = C.S.spk.y;
    
    poshd = C.S.pos.head_direction;
    spkhd = C.S.spk.head_direction;
    
    post = C.S.pos.t;
    spkt = C.S.spk.t;
    
    %% border score
    [border, shuffle_border] = ANA_Border_Score(C, 1, 1/3);
    histogram(shuffle_border, 100);
    hold on;
    line([border, border], ylim, 'LineWidth', 1, 'Color', 'r');
    xlabel('Border Score');
    ylabel('Count');
    
    nless = sum(shuffle_border < border);
    nequal = sum(shuffle_border == border);
    centile = 100 * (nless + 0.5*nequal) / length(shuffle_border);
    
    title(sprintf('Centrile = %.2f', centile));
    
    figfile = fullfile(export_folder, sprintf('%i_ShuffleBorder', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    %% load spike file
    
    Filename = fullfile(C.P.path_dataout, C.P.datadir_out, C.S.spike_file);
    FieldSelectionFlags = [1 1 1 1 1];
    HeaderExtractionFlag = 1;
    ExtractMode = 1;
    ExtractionModeVector = [];
    
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
        Nlx2MatSpike( Filename, FieldSelectionFlags,...
        HeaderExtractionFlag, ExtractMode, ExtractionModeVector);
    
    %% spike shape
    ADBits2uV = 0.00305175781;
    
    for ichan = 1:4
        plot(squeeze(Samples(:, ichan, CellNumbers == C.cell_id).*ADBits2uV), 'Color', [220 220 220]./255, 'LineWidth', .5);
        hold on
        plot(mean(squeeze(Samples(:, ichan, CellNumbers == C.cell_id).*ADBits2uV), 2), 'k', 'LineWidth', 1);
        ylim([-60 100]);
        xlabel('Time (ms)');
        ylabel('Voltage (uV)');
        box off
        
        figfile = fullfile(export_folder, sprintf('%i_SpikeShapes_c%i', cells(ii), ichan));
        print(figfile, '-depsc', '-painters');
        print(figfile, '-dsvg', '-painters');
        print(figfile, '-dpdf', '-painters');
        savefig(figfile);
        close all;
    end
    
    %% behavior
    plot(posx, posy, 'k', spkx, spky, '.r');
    axis tight
    ax = gca;
    ax.Visible = 'off';
    
    figfile = fullfile(export_folder, sprintf('%i_Path', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    %% heatmap
    imsc(C.S.rate_mat, 'jet');
    colorbar
    set(gca, 'Ydir', 'Normal');
    ax = gca;
    ax.Visible = 'off';
    
    figfile = fullfile(export_folder, sprintf('%i_HeatMap', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    %% ISI
    ISI = diff(C.S.spk.t);
    y = ISI;
    x = logspace(-4, 1, 100);
    histogram(y, x);
    set(gca, 'xscale', 'log');
    hold on
    line([0.002 0.002], [0 max(histcounts(log(ISI), 100))], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
    xlabel('Time (ms)');
    ylabel('ISI Count');
    box off
    
    figfile = fullfile(export_folder, sprintf('%i_ISI', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    %% HD/Time
    c_poshd = poshd;
    c_poshd(abs(diff(poshd)) > 300) = NaN;
    
    plot(post - min(post), c_poshd, '-k', spkt - min(post), spkhd, '.r');
    xlabel('Time (sec)');
    ylabel('Direction (degrees)');
    box off
    
    figfile = fullfile(export_folder, sprintf('%i_HeadDirectionxTime', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    % circular histogram
    plot_ax = linspace(-pi, pi, 61);
    
    ax = linspace(-pi, pi, 60);
    y_vec = 0 : C.S.HD.ray_score / 20 : C.S.HD.ray_score; % length of rayleigh vector
    x_vec = ones(1, length(y_vec)) * C.S.HD.ray_angle; % direction of rayleigh vector
    
    polarhistogram('BinEdges', plot_ax, 'BinCounts', C.S.HD.smooth_rate./max(C.S.HD.smooth_rate));
    hold on
    polarplot(x_vec, y_vec, 'r', 'LineWidth', 3);
    hold on
    polarplot(linspace(-pi, pi, 60), smooth(C.S.HD.time_phi)./max(smooth(C.S.HD.time_phi)), '--', 'LineWidth', 3);
    legend('off');
    
    figfile = fullfile(export_folder, sprintf('%i_PolarHD', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    % Rayleigh Shuffle
    
    spk_t_shuffle = ANA_timestamps_shuffle(C, 1000);
    
    HD_rayleigh = zeros(1000, 1);
    
    for ii_idx = 1:1000
        
        shuff_s = C.S;
        shuff_s.spk.t = spk_t_shuffle{ii_idx}';
        shuff_s.spk.head_direction = interp1(shuff_s.pos.t,...
            shuff_s.pos.head_direction,...
            shuff_s.spk.t);
        
        HD_temp = ANA_HD(shuff_s, parms);
        HD_rayleigh(ii_idx) = HD_temp.ray_score;
    end
    
    histogram(HD_rayleigh, 100);
    hold on;
    line([C.S.HD.ray_score, C.S.HD.ray_score], ylim, 'LineWidth', 1, 'Color', 'r');
    xlabel('Rayleigh Score');
    ylabel('Count');
    nless = sum(HD_rayleigh < C.S.HD.ray_score);
    nequal = sum(HD_rayleigh == C.S.HD.ray_score);
    centile = 100 * (nless + 0.5*nequal) / length(HD_rayleigh);
    
    figfile = fullfile(export_folder, sprintf('%i_ShuffleHD', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
    
    %% autocorrelation
    k = dsearchn(C.S.pos.t, C.S.spk.t');
    
    train = zeros(1, length(C.S.pos.t));
    train(k) = 1;
    
    [acor, lags] = xcorr(train, parms.max_lag, 'coeff');
    stem(lags, acor)
    ax = -parms.max_lag : 1 : parms.max_lag;
    
    plot(ax, acor, 'k');
    hold on;
    
    line([-125 -125],[0 max(acor)*0.6],'Color','r')% theta
    line([125 125],[0 max(acor)*0.6],'Color','r')% theta
    line([-250 -250],[0 max(acor)*0.3],'Color','r')%,'LineStyle',':');
    line([250 250],[0 max(acor)*0.3],'Color','r')%,'LineStyle',':');
    line([0 0],[0 max(acor)],'Color','r')%,'LineStyle',':');
    
    xlim([-parms.max_lag parms.max_lag]);
    set(gca,'XTick',[-500 -125 0 125 500]);
    box off;
    xlabel('Lags');
    ylabel('Correlation');
    
    figfile = fullfile(export_folder, sprintf('%i_Autocorrelation', cells(ii)));
    print(figfile, '-depsc', '-painters');
    print(figfile, '-dsvg', '-painters');
    print(figfile, '-dpdf', '-painters');
    savefig(figfile);
    close all;
end

