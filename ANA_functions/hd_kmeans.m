% function hd_kmeans(c, vt, p, s)
clearvars hdScore rLength hdOccupancy hdSpikes hdRates
nHdBins = 60;
numWindows = 3;
windows = linspace(vt.timestamps(1), vt.timestamps(end), numWindows);
for iiWindow = 1:(numWindows-1)
    idxPos = find(vt.timestamps >= windows(iiWindow) & vt.timestamps <= windows(iiWindow+1));
    idxSpk = find(c.timestamps >= windows(iiWindow) & c.timestamps <= windows(iiWindow+1));
    
    windowHdPos = vt.poshd(idxPos);
    windowHdSpk = c.poshd(idxSpk);
    
    hdBins = -pi : 2*pi/nHdBins : pi;
    
    hdOccupancy{iiWindow} = histcounts(windowHdPos, hdBins);
    hdSpikes{iiWindow} = histcounts((windowHdSpk), hdBins);
    hdRates{iiWindow} = hdSpikes{iiWindow} ./ (hdOccupancy{iiWindow} .* (0.04 * 1e-6));
    
    % convert coordinates from polar to rectangular form
    hdBins = -pi : 2*pi/nHdBins : pi - pi/nHdBins;
    x = hdRates{iiWindow} .* cos(hdBins);
    y = hdRates{iiWindow} .* sin(hdBins);
    sumLength = sum(sqrt(x.^2 + y.^2));
    
    % computes sum vector
    rx = sum(x);
    ry = sum(y);
    rLength(iiWindow) = sqrt(rx^2 + ry^2);
    
    hdScore(iiWindow) = rLength(iiWindow) / sumLength;
end
% end