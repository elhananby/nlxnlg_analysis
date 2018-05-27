clearvars;
dbstop if error

idx = 1;
set(groot,...
    'DefaultAxesFontSize', 18,...
    'DefaultAxesFontWeight', 'normal',...
    'DefaultAxesFontName', 'Gill Sans MT',...
    'DefaultPolarAxesFontSize', 18,...
    'DefaultPolarAxesFontWeight', 'normal',...
    'DefaultPolarAxesFontName', 'Gill Sans MT');

[num,txt,raw] = xlsread('elhi_inclusion_list_sessions.xlsx', 'Cells');
cells = txt(2:end,1);

hd_cells_file = fopen('hd_cells.txt', 'a');
proc_cells_file = fopen('proc_cells.txt', 'a');

hd_cells = fileread('hd_cells.txt');
proc_cells = fileread('proc_cells.txt');

% hd_cells = dlmread('hd_cells.txt');
% hd_cells = hd_cells(:,1);

for i = 1:length(cells)
    
    file_search = sprintf('D:\\experiment_data\\cells\\%s_*.mat', num2str(cells{i}));
    file_to_load = subdir(file_search);
    
    if isempty(file_to_load) || contains(proc_cells, num2str(cells{i}))
        continue;
    end
    

    current_cell = file_to_load.name;
    
    % load current cell file
    load(current_cell);
    
    fprintf('Cell #%i\\%i (%.2f%%)\n', i, length(cells), i*100/length(cells));
    
    trial_time = ((C.S.end_time-C.S.start_time)*10^-6)/60;       
    % run HD analysis
    parms.num_of_direction_bins = 60;
    
    HD = ANA_HD(C.S, parms);
    
    if ~exist('C.S.HD', 'var')
        C.S.HD = HD;
    end
    
    %shuffling
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
    
    nless = sum(HD_rayleigh < C.S.HD.ray_score);
    nequal = sum(HD_rayleigh == C.S.HD.ray_score);
    centile = 100 * (nless + 0.5*nequal) / length(HD_rayleigh);
    
    fprintf(proc_cells_file, '%s\r\n', num2str(cells{i}));
    
    if ~(trial_time >= 10 && centile >= 95 && length(C.S.spk.t) >= 300 && C.S.HD.ray_pval <= 0.1)
        continue;
    end
    
    fprintf('\tTime = %.2f\tCentile = %.2f\tpval = %.2f\n', trial_time, centile, C.S.HD.ray_pval);
    
    % setup for half analysis
    start_t(1) = C.S.pos.t(1); end_t(1) = C.S.pos.t(floor(end/2));
    start_t(2) = C.S.pos.t(floor(end/2) + 1); end_t(2) = C.S.pos.t(end);
    
    rate_ang_smooth(idx, :) = HD.smooth_rate; % get smoothed spike rate histogram
    
    if length(num2str(C.day)) == 6
        date = datetime(num2str(C.day), 'InputFormat', 'ddMMyy');
    elseif length(num2str(C.day)) == 5
        date = datetime(num2str(C.day), 'InputFormat', 'dMMyy');
    end
    
    if date.Year == 2017 % if 2017, use old scripts
        rate_ang_smooth(idx, :) = circshift(rate_ang_smooth(idx, :), 15);
    end
    
    for ii_half = 1:2
        time_idx = C.S.pos.t >= start_t(ii_half) & C.S.pos.t <= end_t(ii_half);
        spk_idx = C.S.spk.t >= start_t(ii_half) & C.S.spk.t <= end_t(ii_half);
        
        dt = mean(diff(C.S.pos.t));
        n_bins = 60;
        
        HD_time_phi = deg2rad(wrapTo180(C.S.pos.head_direction(time_idx)));
        HD_count_phi = deg2rad(wrapTo180(C.S.spk.head_direction(spk_idx)));
        
        ang_ax = linspace(-pi, pi, 61);
        count_phi = histcounts(HD_count_phi, ang_ax);
        time_phi = histcounts(HD_time_phi, ang_ax) * dt;
        rate_ang = count_phi./time_phi;
        
        ana_ax = linspace(-pi, pi, 60);
        HD_ray_angle(idx, ii_half) = circ_mean(ana_ax', rate_ang');
    end
    
    animal_num = split(file_to_load.folder, '\');
    animal_num = split(animal_num{4}, '_');
    
    animal_list(idx) = str2double(animal_num{1});   
    
    idx = idx + 1;
    fprintf(hd_cells_file, '%s\t%i\t%i\r\n', num2str(cells{i}), C.animal, C.day);
end

fclose('all');

animal_list(any(isnan(HD_ray_angle), 2)) = [];
HD_ray_angle(any(isnan(HD_ray_angle), 2), :) = [];

hd_half(HD_ray_angle, animal_list');
hd_variance(rate_ang_smooth, animal_list);

function hd_half(half_ray_angle, animal_color)

cd('D:\');
[rho, pval] = circ_corrcc(half_ray_angle(:,1), half_ray_angle(:,2)); % circular correlation

hd_angle_deg = wrapTo360(rad2deg(half_ray_angle)); % convert to 360-degrees

hd_angle_fix = hd_angle_deg; % initialize

for ii = 1:length(hd_angle_deg)
    
    linear_diff = abs(hd_angle_deg(ii, 1) - hd_angle_deg(ii, 2)); % find difference between two points
    
    if linear_diff > 180 % if linear difference is bigger than 180 degrees
        
        if hd_angle_deg(ii, 1) > hd_angle_deg(ii, 2) % if the first half degree is larger than second half
            
            angle_diff = (360 - hd_angle_deg(ii, 1)) + hd_angle_deg(ii, 2); % find the angular difference
            hd_angle_fix(ii, 2) = hd_angle_deg(ii, 1) + angle_diff; % add it to the lower value
        else
            angle_diff = (360 - hd_angle_deg(ii, 2)) + hd_angle_deg(ii, 1);
            hd_angle_fix(ii, 1) = hd_angle_deg(ii, 2) + angle_diff;
        end
        
    end
    
end

figure;
gscatter(hd_angle_fix(:, 1), hd_angle_fix(:, 2), animal_color, [], [], 20);
corr_line = polyfit(hd_angle_fix(:, 1), hd_angle_fix(:, 2), 1);
y_pred = polyval(corr_line, hd_angle_fix(:,1));
hold on;
plot(hd_angle_fix(:, 1), y_pred, 'Color', 'black', 'LineWidth', 1);

xlabel('1st Half Rayleigh Angle', 'FontSize', 24, 'FontName', 'Gill Sans MT');
ylabel('2st Half Rayleigh Angle', 'FontSize', 24, 'FontName', 'Gill Sans MT');
%title(sprintf('rho = %.2f, p = %.2e', rho, pval), 'FontSize', 24);
legend('off');
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20.9 15.5191],...
    'Color', 'w');
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20.9 15.5191]);
set(gca,'LooseInset',get(gca,'TightInset')) 
export_fig(sprintf('hdcorr'), '-png', '-r300', '-cmyk', '-transparent', '-nofontswap');
close all;
end

function hd_variance(rate_ang_smooth, animal_list)
cd('D:\');
[m, n] = size(rate_ang_smooth);
rate_ang_norm = zeros(m,n);

for ii = 1:m
    rate_ang_norm(ii, :) = rescale(rate_ang_smooth(ii, :));
end

keep_idx = rate_ang_norm > 0.8;
[sort_animal, order_animal] = sort(animal_list);
new_keep_idx = keep_idx(order_animal, :);
colored_hd = bsxfun(@times, new_keep_idx, sort_animal');

%% singular heatmap
subplot(10,2,[1,18])
imagesc(colored_hd);
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca, 'box', 'off',...
    'XTickLabel', [], 'XTick', [], 'XColor', 'w',...
    'FontSize', 18, 'FontName', 'Gill Sans MT');
ylabel('Cell Number', 'FontSize', 24, 'FontName', 'Gill Sans MT');


subplot(10,2,[19,20])
histogram('BinEdges', 0:1:60, 'BinCounts', sum(keep_idx, 1));
set(gca,'LooseInset',get(gca,'TightInset')) 

xlim([0 60]);
xticks([0 15 30 45 60]);
xticklabels({'-180', '-90', '0', '90', '180'});
set(gca,'LooseInset',get(gca,'TightInset')) 
xlabel('Head-direction (degrees)', 'FontSize', 24, 'FontName', 'Gill Sans MT');
set(gca, 'box', 'off',...
    'YTick', [], 'YTickLabel', [], 'YColor', 'w')
set(gca,'LooseInset',get(gca,'TightInset')) 
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20.9 16],...
    'Color', 'w');

set(gca,'LooseInset',get(gca,'TightInset')) 
export_fig(sprintf('hdvariance'), '-png', '-r300', '-cmyk');

hold off;

close all

end