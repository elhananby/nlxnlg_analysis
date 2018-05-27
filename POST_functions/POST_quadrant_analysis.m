function ANA_quadrant_analysis( c )
%ANA_QUADRANT_ANALYSIS Run all basic analysis by arena quadrant
%   INPUT:  c - cell output .mat file
%   OUTPUT: Figure per quadrant

%parameters used for analysis
parms.sigma = 3; % gaussiam smoothing factor
parms.time_per_bin= 0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
parms.bin_size = 3; % size of spacial bin (for create the rate map)
parms.num_of_direction_bins = 60; % for head-direction calculation
parms.max_lag = 500; % max lag (in msec) for temporal autocorrelation

q = c.S; % create new variable q for
min_x = min([q.VT_ArenaMarking.x_NW,q.VT_ArenaMarking.x_NE,q.VT_ArenaMarking.x_SE,q.VT_ArenaMarking.x_SW]);
max_x = max([q.VT_ArenaMarking.x_NW,q.VT_ArenaMarking.x_NE,q.VT_ArenaMarking.x_SE,q.VT_ArenaMarking.x_SW]);
min_y = min([q.VT_ArenaMarking.y_NW,q.VT_ArenaMarking.y_NE,q.VT_ArenaMarking.y_SE,q.VT_ArenaMarking.y_SW]);
max_y = max([q.VT_ArenaMarking.y_NW,q.VT_ArenaMarking.y_NE,q.VT_ArenaMarking.y_SE,q.VT_ArenaMarking.y_SW]);

for ii_quadrant = 1:4
    q = c.S;
    
    if ii_quadrant == 1 % NW
        q.VT_ArenaMarking.min_y = min_y + (max_y - min_y)/2;
        q.VT_ArenaMarking.max_y = max_y;
        q.VT_ArenaMarking.min_x = min_x;
        q.VT_ArenaMarking.max_x = max_x - (max_x - min_x)/2;
    elseif ii_quadrant == 2 % NE
        q.VT_ArenaMarking.min_y = min_y + (max_y - min_y)/2;
        q.VT_ArenaMarking.max_y = max_y;
        q.VT_ArenaMarking.min_x = min_x + (max_x - min_x)/2;
        q.VT_ArenaMarking.max_x = max_x;
    elseif ii_quadrant == 3 % SE
        q.VT_ArenaMarking.min_y = min_y;
        q.VT_ArenaMarking.max_y = max_y - (max_y - min_y)/2;
        q.VT_ArenaMarking.min_x = min_x + (max_x - min_x)/2;
        q.VT_ArenaMarking.max_x = max_x;
    elseif ii_quadrant == 4 % SW
        q.VT_ArenaMarking.min_y = min_y;
        q.VT_ArenaMarking.max_y = max_y - (max_y - min_y)/2;
        q.VT_ArenaMarking.min_x = min_x;
        q.VT_ArenaMarking.max_x = max_x - (max_x - min_x)/2;
    end
    
    q_pos_idx = find(q.pos.x >= q.VT_ArenaMarking.min_x &...
        q.pos.x <= q.VT_ArenaMarking.max_x &...
        q.pos.y >= q.VT_ArenaMarking.min_y &...
        q.pos.y <= q.VT_ArenaMarking.max_y);
    
    q_spk_idx = find(q.spk.x >= q.VT_ArenaMarking.min_x &...
        q.spk.x <= q.VT_ArenaMarking.max_x &...
        q.spk.y >= q.VT_ArenaMarking.min_y &...
        q.spk.y <= q.VT_ArenaMarking.max_y);
    
    % apply indexing to pos data
    q.pos.t = q.pos.t(q_pos_idx);
    q.pos.x = q.pos.x(q_pos_idx);
    q.pos.y = q.pos.y(q_pos_idx);
    q.pos.x1 = q.pos.x1(q_pos_idx);
    q.pos.y1 = q.pos.y1(q_pos_idx);
    q.pos.x2 = q.pos.x2(q_pos_idx);
    q.pos.y2 = q.pos.y2(q_pos_idx);
    q.pos.head_direction = q.pos.head_direction(q_pos_idx);
    q.pos.vx = q.pos.vx(q_pos_idx);
    q.pos.vy = q.pos.vy(q_pos_idx);
    
    % apply indexing to spk data
    q.spk.t = q.spk.t(q_spk_idx);
    q.spk.x = q.spk.x(q_spk_idx);
    q.spk.y = q.spk.y(q_spk_idx);
    q.spk.x1 = q.spk.x1(q_spk_idx);
    q.spk.y1 = q.spk.y1(q_spk_idx);
    q.spk.x2 = q.spk.x2(q_spk_idx);
    q.spk.y2 = q.spk.y2(q_spk_idx);
    q.spk.head_direction = q.spk.head_direction(q_spk_idx);
    q.spk.vx = q.spk.vx(q_spk_idx);
    q.spk.vy = q.spk.vy(q_spk_idx);
    
    q = ANA_Common_Spatial_Calc(q, parms, 1);
    q.q_flag = 1;
    q.quadrant = ii_quadrant;
    fig = ANA_Plot_figure(c, q, parms, c.cell_number, 1);
    
    file_name = sprintf('%s_q%d', c.cell_file_name(1:end-4), ii_quadrant);
    save_file_name = fullfile(c.P.path_dataout, 'cells', file_name);
    saveas(fig, save_file_name, 'jpg'); 
    %savefig(fig_1, save_file_name, 'compact');
    fprintf('\t\tsaved figure %s\n',save_file_name);
end
end

