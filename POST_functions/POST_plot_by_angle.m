function POST_plot_by_angle(c)

[~, idx] = max(movmean(c.HD.rate_phi, n_bins/2));
ax = -180 : 360 / (n_bins - 1) : 180;
val = ax(idx);
ana_idx = cell(2,1);
if (ax(idx) <= 90 && ax(idx) >= -90)
    range = [ax(idx - 15) ax(idx + 15)];
    ana_idx{1} = c.pos.head_direction > range(1) & c.pos.head_direction < range(2);
    ana_idx{2} = not(ana_idx{1});
elseif ax(idx) < -90
    inv_range =  [ax(idx + 15) ax(idx + 45)];
    ana_idx{2} = c.pos.head_direction > inv_range(1) & c.pos.head_direction < inv_range(2);
    ana_idx{1} = not(ana_idx{2});
elseif ax(idx) > 90
    inv_range = [ax(idx - 15) ax(idx - 45)];
    ana_idx{2} = c.pos.head_direction > inv_range(1) & c.pos.head_direction < inv_range(2);
    ana_idx{1} = not(ana_idx{2});
end

for ii = 1:2
    angle_S = ANA_Common_Spatial_Calc(S_range, parms, 0);
    fig = ANA_Plot_figure(C, angle_S, parms, c.cell_number,session_num);
    file_name = sprintf('maxangle%d_%s_%.2d_%s', ii, C.cell_file_name(1:end-4), session_num, c.session);
    save_file_name = fullfile(C.P.path_dataout, 'cells', file_name);
end

end
