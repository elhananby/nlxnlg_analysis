clearvars
close all
dbstop if error

cells = importdata('hdcells.txt');

% col 1 is file name, col 2 is hdOccupancy, col 3 is hdRates, col 4 is hdScore
% col 5 is mean direction, col 6 hdscore1, col 7 is hdscore2
cellData = cell(length(cells), 7);

idx = 1;
animal_list = [];


for ii = 1:numel(cells)
        
    % load file
    files_search = sprintf('D:\\experiment_data\\cells\\%s_*.mat', num2str(cells(ii)));
    file_to_load = subdir(files_search);
    load(file_to_load.name);
    
    [~, name, ~] = fileparts(file_to_load.name);
    number = regexp(name, 'quail(\d+)', 'tokens');
    animal_list(ii) = str2double(cell2mat(number{1}));
    
    disp(['Cell ', num2str(cells(ii)), ' (',  num2str(ii*100/numel(cells)), '%)']);
    
    % calculate hd data
    t = C.S.pos.t;
    ts = C.S.spk.t';
    dt = mean(diff(t));
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)';
    direction = wrapTo2Pi(deg2rad(C.S.pos.head_direction));
    
    % save cell data
    cellData{ii, 1} = file_to_load.name;
    [cellData{ii, 2}, cellData{ii, 3}, cellData{ii, 4}, cellData{ii, 5}] = calculate_hd_score(direction, spiketrain, dt, 60);
    
   
    half_idx = ceil(length(t)/2);
    half_time = t(half_idx);
    
    % first half    
    t = C.S.pos.t(C.S.pos.t <= half_time );
    ts = C.S.spk.t(C.S.spk.t <= half_time)';
    
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)';
    direction = wrapTo2Pi(deg2rad(C.S.pos.head_direction(C.S.pos.t <= half_time)));
    
    [~, ~, ~, cellData{ii, 6}] = calculate_hd_score(direction, spiketrain, dt, 60);
    
    % second half
    t = C.S.pos.t(C.S.pos.t > half_time);
    ts = C.S.spk.t(C.S.spk.t > half_time)';
    
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)';
    direction = wrapTo2Pi(deg2rad(C.S.pos.head_direction(C.S.pos.t > half_time)));
    
    [~, ~, ~, cellData{ii, 7}] = calculate_hd_score(direction, spiketrain, dt, 60);
end
%%
hd_ratemap(cellData)
% hd_correlation(cellData, animal_list)

function hd_ratemap(cellData)
export_folder = 'D:\Dropbox (Technion Dropbox)\Lab_folder\Grant_2018_data';
meanHD = cell2mat(cellData(:,5));
meanHD = wrapTo2Pi(meanHD);
[~, I] = sort(meanHD); % sort by mean direction

Win = hamming(10);
Win = Win/sum(Win);
rate_mat = nan(length(cellData), 60);

for ii = 1:length(cellData)
    temp_mat = normalize(cellData{ii, 3}');
    rate_mat(ii, :) = cconv(temp_mat,...
        Win, 60);
end

imagesc(rate_mat(I, :));
colormap(jet);
xlabel('Direction');
ylabel('Cell #');
xlim([0 60]);
xticks([1 15 30 45 60]);
xticklabels({'-180', '-90', '0', '90', '180'});
box off
axis tight

figfile = fullfile(export_folder, sprintf('hd_population'));
print(figfile, '-depsc', '-painters');
print(figfile, '-dsvg', '-painters');
print(figfile, '-dpdf', '-painters');
savefig(figfile);
close all;

end
function hd_correlation(cellData, animal_list)

export_folder = 'D:\Dropbox (Technion Dropbox)\Lab_folder\Grant_2018_data';
hd_half = [cell2mat(cellData(:, 6)), cell2mat(cellData(:, 7))];
[rho, pval] = circ_corrcc( cell2mat(cellData(:, 6)), cell2mat(cellData(:, 7))  );

hd_angle_deg = wrapTo360(rad2deg(hd_half)); % convert to 360-degrees

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
gscatter(hd_angle_fix(:, 1), hd_angle_fix(:, 2), animal_list, brewermap(unique(length(animal_list)), 'Paired'));

hold on;
xl = xlim;
yl = ylim;
line(xl, yl, 'Color', 'black')
xlabel('Mean Direction (1st Half)');
ylabel('Mean Direction (2nd Half)');
legend off

figfile = fullfile(export_folder, sprintf('hd_correlation'));
print(figfile, '-depsc', '-painters');
print(figfile, '-dsvg', '-painters');
print(figfile, '-dpdf', '-painters');
savefig(figfile);
close all;

end