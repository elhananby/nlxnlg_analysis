
dbstop if error

clearvars;
close all;

% fid = fopen('cells.txt', 'rt');

s = tdfread('D:\hdcells.txt');
cells = s.cell_num(logical(s.use));

%% poster plotting
fig_size = 6; % in cm
set(groot,...
    'DefaultAxesFontSize', 16,...
    'DefaultAxesFontWeight', 'normal',...
    'DefaultAxesFontName', 'Gill Sans MT',...
    'DefaultPolarAxesFontSize', 16,...
    'DefaultPolarAxesFontWeight', 'normal',...
    'DefaultPolarAxesFontName', 'Gill Sans MT');

cd('D:\OneDrive - Technion\PhD\imprs poster\cells')

%while true
for ii = 1:numel(cells)
%     cell_num = fgetl(fid);
%     
%     if ~ischar(cell_num)
%         disp('EOF');
%         break;
%     end

    fprintf('Processing cell %s...\n', num2str(cells(ii)));
    file_search = sprintf('D:\\experiment_data\\cells\\%s_*.mat', num2str(cells(ii)));
    file_to_load = subdir(file_search);
    
    load(file_to_load.name);
    
    % clean out-of-bounds
    minx = C.S.VT_ArenaMarking.min_x;
    maxx = C.S.VT_ArenaMarking.max_x;
    miny = C.S.VT_ArenaMarking.min_y;
    maxy = C.S.VT_ArenaMarking.max_y;
    
    outx = find(C.S.pos.x >= maxx | C.S.pos.x <= minx);
    outy = find(C.S.pos.y >= maxy | C.S.pos.y <= miny);
    out = [outx; outy];
    out = unique(out);
    C.S.pos.x(out) = NaN;
    C.S.pos.y(out) = NaN;
    
    outx = find(C.S.spk.x >= maxx | C.S.spk.x <= minx);
    outy = find(C.S.spk.y >= maxy | C.S.spk.y <= miny);
    [m,n] = size(outx);
    if n > m
        outx = outx';
        outy = outy';
    end
    
    try
        out = [outx; outy];
    catch
        out = [outx, outy];
    end
        
    out = unique(out);
    C.S.spk.x(out) = NaN;
    C.S.spk.y(out) = NaN;
    
    %% behavior
    figure;
    axes('Units', 'normalized', 'Position', [0 0 1 1]);
    plot(C.S.pos.x - min(C.S.pos.x), C.S.pos.y - min(C.S.pos.y), 'k', C.S.spk.x - min(C.S.pos.x), C.S.spk.y - min(C.S.pos.y), '.r');
    
    axis off
    box off
    
    set(gca,'YDir','reverse')
    
    set(gcf, 'units', 'centimeters'	, 'position', [0 0 fig_size fig_size])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    filename = sprintf('%d_behavior', C.cell_number);
    export_fig(filename, '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap')
    export_fig(filename, '-eps', '-nofontswap')
    savefig(filename)
    
    %% rate map
    %close all;
    figure;
    dt = mean(diff(C.S.pos.t));
    [time, xedges, yedges] = histcounts2(C.S.pos.x, C.S.pos.y, 30);
    time = time .* dt;
    
    count = histcounts2(C.S.spk.x, C.S.spk.y, xedges, yedges);
    rate = count ./ time;
    clean_rate = rate;
    clean_rate(clean_rate == 0) = NaN;
    [tf, down, up, center] = isoutlier(clean_rate(:));
    rate(rate >= up) = NaN;
    rate(isinf(rate)) = NaN;
    rate_mat_nanconv = nanconv(rate, fspecial('gaussian', 3*[3,3], 3), 'edge');
    
    sc(rate_mat_nanconv, jet)
    
    box off
    axis off
    firing_rate = length(C.S.spk.t)/(C.S.pos.t(end) - C.S.pos.t(1));
     
    set(gcf, 'units', 'centimeters'	, 'position', [0 0 fig_size fig_size])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    filename = sprintf('%d_rate', C.cell_number);
    export_fig(filename, '-png', '-r300', '-cmyk', '-nofontswap')
    export_fig(filename, '-eps', '-nofontswap')
    savefig(filename)
    %% head-direction
    
    figure;
    dt = mean(diff(C.S.pos.t));
    ax = linspace(-pi, pi, 60);
    hd_time = histcounts(C.S.pos.head_direction, 60) .* dt;
    hd_count = histcounts(C.S.spk.head_direction, 60);
    hd_rate = hd_count./hd_time;
    
    Win = hamming(10);
    Win = Win/sum(Win);
    smooth_rate = cconv(hd_rate, Win, 60);
    %norm_hd_time = rescale(smooth_time);
    
    ray_score = circ_r(ax', smooth_rate');
    ray_angle = circ_mean(ax', smooth_rate');
    
    %     polarplot(ax, norm_hd_time, '-ob', ax, smooth_rate, '.-r');
    
    plot_ax = linspace(-pi, pi, 61);
    polarhistogram('BinEdges', plot_ax, 'BinCounts', smooth_rate./max(smooth_rate))
    hold on;
    y_vec = 0 : ray_score / 20 : ray_score; % length of rayleigh vector
    x_vec = ones(1, length(y_vec)) * ray_angle; % direction of rayleigh vector
    polarplot(x_vec, y_vec, 'r', 'LineWidth', 2);
    
    
    set(gcf, 'Units', 'Centimeters', 'Position', [0 0 fig_size fig_size]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    filename = sprintf('%d_polar', C.cell_number);
    export_fig(filename, '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap')
    export_fig(filename, '-eps', '-nofontswap')
    savefig(filename)
    hold off;
    
    %% hd/time
    figure;
    pos_time_vector = (C.S.pos.t - min(C.S.pos.t))/60;
    spk_time_vector = (C.S.spk.t - min(C.S.pos.t))/60;
 
    plot(pos_time_vector, C.S.pos.head_direction,...
        spk_time_vector, C.S.spk.head_direction, '.r');
    
    ylim([-180 180]);
    xlim([min(pos_time_vector), max(pos_time_vector)]);
    
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Head Direction', 'FontSize', 18);
    
    set(gcf, 'units', 'centimeters', 'position', [0 0 fig_size*2 fig_size]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    filename = sprintf('%d_hdtime', C.cell_number);
    export_fig(filename, '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap')
    export_fig(filename, '-eps', '-nofontswap')
    savefig(filename)
    %% rayleigh shuffle
    figure;
    parms.num_of_direction_bins = 60;
    
    C.S.HD = ANA_HD(C.S, parms);
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
    
    histogram(HD_rayleigh, 100, 'EdgeColor', 'none');
    hold on;
    line([C.S.HD.ray_score, C.S.HD.ray_score], ylim, 'LineWidth', 1, 'Color', 'r');
    xlabel('Rayleigh Score', 'FontSize', 18);
    ylabel('Count', 'FontSize', 18);

    set(gcf, 'Units', 'centimeters', 'Position', [0 0 fig_size*1.3 fig_size]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    filename = sprintf('%d_rayleigh', C.cell_number);
    export_fig(sprintf('%s', filename), '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap');
    export_fig(sprintf('%s', filename), '-eps', '-nofontswap')
    savefig(sprintf('%s', filename))
    
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 fig_size fig_size]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    filename = sprintf('%d_rayleigh-square', C.cell_number);
    export_fig(sprintf('%s', filename), '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap');
    export_fig(sprintf('%s', filename), '-eps', '-nofontswap')
    savefig(sprintf('%s', filename))
    hold off;
    
    %% ISI
    figure;
    ISI = diff(C.S.spk.t);
    ax = 0 : 0.001 : max(ISI);
    h = histogram(ISI, ax, 'LineStyle', 'none', 'EdgeColor', 'none');
    set(gca, 'xscale', 'log');
    %h.Normalization = 'countdensity';
    hold on
    line([0.002 0.002], [0 max(histcounts(ISI, ax))], 'Color', 'red', 'LineStyle', '--')
    set(findobj(gcf, 'Type', 'line'), 'LineWidth', 2);
    xlim([0 1]);
    
    box off
    
    xlabel('Time (Seconds)', 'FontSize', 18);
    ylabel('Count', 'FontSize', 18)
    
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 fig_size*1.3 fig_size]);
    filename = sprintf('%d_isi', C.cell_number);
    export_fig(sprintf('%s', filename), '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap');
        export_fig(filename, '-eps', '-nofontswap')
        savefig(filename)
    hold off
    
    %% spike shape
    figure;
    time_ax = linspace(1, 1000, 32); % there are 32 time bins
    SpikeShape = (C.S.mean_spike_shape)*0.01; %voltage from NLX is in uv*100
    plot(time_ax, SpikeShape(:,1),'k'); hold on;
    plot(time_ax, SpikeShape(:,2),'r');
    plot(time_ax, SpikeShape(:,3),'g');
    plot(time_ax, SpikeShape(:,4),'b');
    
    box off;
    
    ylabel('Voltage (uV)', 'FontSize', 18);
    xlabel('Time (uS)', 'FontSize', 18)
    
      set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'units', 'centimeters'	, 'position', [0 0 fig_size fig_size]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    filename = sprintf('%d_shape', C.cell_number);
    export_fig(filename, '-png', '-transparent',  '-r300', '-cmyk', '-nofontswap')
        export_fig(filename, '-eps', '-nofontswap')
        savefig(filename)
    
    %%
    close all
end
