function [borderOccupancy, borderSpikes, borderRates] = calculate_border_map(c, vt)
global boxSize dt;
distVec = 0:boxSize/10:boxSize;

%% X axis (east-west)
borderOccupancy{1} = histcounts(vt.posx_c, distVec);
borderSpikes{1} = histcounts(c.posx_c, distVec);

borderRates{1} = borderSpikes{1} ./ (borderOccupancy{1}.*dt);

%% Y axis (south-north)
borderOccupancy{2} = histcounts(vt.posy_c, distVec);
borderSpikes{2} = histcounts(c.posy_c, distVec);

borderRates{2} = borderSpikes{2} ./ (borderOccupancy{2}.*dt);

end

