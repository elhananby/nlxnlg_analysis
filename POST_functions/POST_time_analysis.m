function OUT_time_analysis( c )
%ANA_TIME_ANALYSIS Cut experiment in half and analyse both halves
%   Detailed explanation goes here
%parameters used for analysis
parms.sigma = 3; % gaussiam smoothing factor
parms.time_per_bin= 0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
parms.bin_size = 3; % size of spacial bin (for create the rate map)
parms.num_of_direction_bins = 60; % for head-direction calculation
parms.max_lag = 500; % max lag (in msec) for temporal autocorrelation

start_t = [];
end_t = [];
start_t(1) = c.S.pos.Timestamps(1); end_t(1) = c.S.pos.Timestamps(floor(end/2));
start_t(2) = c.S.pos.Timestamps(floor(end/2) + 1); end_t(2) = c.S.pos.Timestamps(end);



for ii_time = 1:2
    t = c.S;
    t.VT_ArenaMarking.min_x = min([t.VT_ArenaMarking.x_NW,t.VT_ArenaMarking.x_NE,t.VT_ArenaMarking.x_SE,t.VT_ArenaMarking.x_SW]);
    t.VT_ArenaMarking.max_x = max([t.VT_ArenaMarking.x_NW,t.VT_ArenaMarking.x_NE,t.VT_ArenaMarking.x_SE,t.VT_ArenaMarking.x_SW]);
    t.VT_ArenaMarking.min_y = min([t.VT_ArenaMarking.y_NW,t.VT_ArenaMarking.y_NE,t.VT_ArenaMarking.y_SE,t.VT_ArenaMarking.y_SW]);
    t.VT_ArenaMarking.max_y = max([t.VT_ArenaMarking.y_NW,t.VT_ArenaMarking.y_NE,t.VT_ArenaMarking.y_SE,t.VT_ArenaMarking.y_SW]);
    t.q_flag = 0;
    
    t_pos_idx = find(t.pos.Timestamps >= start_t(ii_time) & t.pos.Timestamps <= end_t(ii_time));
    t_spk_idx = find(t.spk.Timestamps >= start_t(ii_time) & t.spk.Timestamps <= end_t(ii_time)); 
   
    % apply indexing to pos data
    t.pos.t = t.pos.t(t_pos_idx);
    t.pos.x = t.pos.x(t_pos_idx);
    t.pos.y = t.pos.y(t_pos_idx);
    t.pos.x1 = t.pos.x1(t_pos_idx);
    t.pos.y1 = t.pos.y1(t_pos_idx);
    t.pos.x2 = t.pos.x2(t_pos_idx);
    t.pos.y2 = t.pos.y2(t_pos_idx);
    t.pos.head_direction = t.pos.head_direction(t_pos_idx);
    t.pos.vx = t.pos.vx(t_pos_idx);
    t.pos.vy = t.pos.vy(t_pos_idx);
    
    % apply indexing to spk data
    t.spk.t = t.spk.t(t_spk_idx);
    t.spk.x = t.spk.x(t_spk_idx);
    t.spk.y = t.spk.y(t_spk_idx);
    t.spk.x1 = t.spk.x1(t_spk_idx);
    t.spk.y1 = t.spk.y1(t_spk_idx);
    t.spk.x2 = t.spk.x2(t_spk_idx);
    t.spk.y2 = t.spk.y2(t_spk_idx);
    t.spk.head_direction = t.spk.head_direction(t_spk_idx);
    t.spk.vx = t.spk.vx(t_spk_idx);
    t.spk.vy = t.spk.vy(t_spk_idx);
    
    t = ANA_Common_Spatial_Calc(t, parms, 1);
    
    fig = ANA_Plot_figure(c, t, parms, c.cell_number, 1);
    file_name = sprintf('%s_t%d', c.cell_file_name(1:end-4), ii_time);
    save_file_name = fullfile(c.P.path_dataout, 'cells', file_name);
    saveas(fig, save_file_name, 'jpg'); 
    fprintf('\t\tsaved figure %s\n',save_file_name);
end
        
end

