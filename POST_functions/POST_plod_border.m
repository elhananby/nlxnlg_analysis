% find threshold values for statistical analysis
borderThr = quantile(shuffledBorderScore, 0.99);
borderThrUnclean = quantile(shuffledBorderScoreUnclean, 0.99);

if orgBorderScore >= borderThr || orgBorderScoreUnclean >= borderThrUnclean
    
    % save to array
    borderCells(borderIdx).cell = c.cell_number;
    borderCells(borderIdx).animal = c.animal;
    borderCells(borderIdx).score = orgBorderScore;
    borderCells(borderIdx).dir = orgDir;
    borderCells(borderIdx).heatmap = orgPosRatesSmooth;
    borderCells(borderIdx).heatmapUnclean = orgPosRatesSmoothUnclean;
    borderIdx = borderIdx + 1;
    
    % save figure
    borderFig = figure('Units', 'pixels', 'Position', [0 0 600 900]);
    
    % behavior
    subplot(5,2,[1 2])
    plot(vt.posx_c, vt.posy_c, '-k', c.posx_c, c.posy_c, '.r');
    switch orgDir
        case 1, title('Border west');
        case 2, title('Border south');
        case 3, title('Border east');
        case 4, title('Border north');
    end
    axis tight
    axis off
    
    %% clean
    % behavior map (cleaned)
    subplot(5,2,3)
    imagesc(orgPosOccupancy);
    colormap(jet);
    title('cleaned (mean + 6std)');
    axis off
    set(gca, 'XDir', 'reverse');
    
    % spike map (cleaned)
    subplot(5,2,5)
    imagesc(orgPosSpikes);
    colormap(jet);
    axis off
    set(gca, 'XDir', 'reverse');
       
    % rate map (cleaned)
    subplot(5,2,7)
    imagesc(orgPosRatesSmooth);
    colormap(jet);
    axis off
    set(gca, 'XDir', 'reverse');
       
    %% original
    % behavior map (cleaned)
    subplot(5,2,4)
    imagesc(orgPosOccupancy);
    colormap(jet);
    title('origianl');
    axis off
    set(gca, 'XDir', 'reverse');
       
    % spike map (cleaned)
    subplot(5,2,6)
    imagesc(orgPosSpikes);
    colormap(jet);
    axis off
    set(gca, 'XDir', 'reverse');
    
    % rate map (cleaned)
    subplot(5,2,8)
    imagesc(orgPosRatesSmoothUnclean);
    colormap(jet);
    axis off
    set(gca, 'XDir', 'reverse');
    
    % shuffling
    subplot(5,2,9)
    histogram(shuffledBorderScore);
    hold on
    line([orgBorderScore orgBorderScore], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    
    subplot(5,2,10)
    histogram(shuffledBorderScoreUnclean);
    hold on
    line([orgBorderScoreUnclean orgBorderScoreUnclean], ylim, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-');
    
    % save figure
    borderFigFile = fullfile(figDir, sprintf('%i_%i_%i_border', c.cell_number, c.animal, c.day));
    export_fig(borderFig, borderFigFile, '-jpg');
    close all;
end