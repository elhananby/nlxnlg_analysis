function PLOT_by_angle(c)
global m n count; % define global variables for plotting
m = 2;
n = 3;
count = 1;

parms.sigma = 3; % gaussiam smoothing factor
parms.time_per_bin= 0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
parms.bin_size = 3; % size of spacial bin (for create the rate map)
parms.num_of_direction_bins = 60; % for head-direction calculation
parms.max_lag = 500; % max lag (in msec) for temporal autocorrelation

S = ANA_Common_Spatial_Calc(c.S, parms, 0);

% calculate a 180-degree range moving mean over each point of rate_phi
[~, idx] = max(movmean(S.HD.rate_phi, parms.num_of_direction_bins/2));
ax = -180 : 360 / (parms.num_of_direction_bins - 1) : 180;
val = ax(idx);
ana_idx = cell(2,1);

[pks, locs] = findpeaks(S.HD.rate_phi);
if (ax(idx) <= 90 && ax(idx) >= -90)
    range = [ax(idx - 15) ax(idx + 15)];
    ana_idx{1} = S.pos.head_direction > range(1) & S.pos.head_direction < range(2);
    ana_idx{2} = not(ana_idx{1});
elseif ax(idx) < -90
    inv_range =  [ax(idx + 15) ax(idx + 45)];
    ana_idx{2} = S.pos.head_direction > inv_range(1) & S.pos.head_direction < inv_range(2);
    ana_idx{1} = not(ana_idx{2});
elseif ax(idx) > 90
    inv_range = [ax(idx - 15) ax(idx - 45)];
    ana_idx{2} = S.pos.head_direction > inv_range(1) & S.pos.head_direction < inv_range(2);
    ana_idx{1} = not(ana_idx{2});
end

for ii = 1:2
    S_range = S;
    S_range.pos.Timestamps = S_range.pos.Timestamps(ana_idx{ii});
    S_range.pos.t = S_range.pos.t(ana_idx{ii});
    S_range.pos.x1 = S_range.pos.x1(ana_idx{ii});
    S_range.pos.x2 = S_range.pos.x2(ana_idx{ii});
    S_range.pos.y1 = S_range.pos.y1(ana_idx{ii});
    S_range.pos.y2 = S_range.pos.y2(ana_idx{ii});
    S_range.pos.x = S_range.pos.x(ana_idx{ii});
    S_range.pos.y = S_range.pos.y(ana_idx{ii});
    S_range.pos.head_direction = S_range.pos.head_direction(ana_idx{ii});
    S_range.pos.vx = S_range.pos.vx(ana_idx{ii});
    S_range.pos.vy = S_range.pos.vy(ana_idx{ii});
    S_range.q_flag = 0;
    
    S_range_new = ANA_Common_Spatial_Calc(S_range, parms, 0);
    
    subplot(m,n,count);
    histogram(S_range_new.pos.head_direction);
    count = count + 1;
    
    PLOT_behavior_map(S_range_new, parms);
    PLOT_rate_map(S_range_new, c);
end

end