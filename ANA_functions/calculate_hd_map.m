function [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt, cKeepIdx, vtKeepIdx)
global nHdBins dt;

hdBins = -pi : 2*pi/nHdBins : pi;

hdOccupancy = histcounts(vt.poshd(vtKeepIdx), hdBins);
hdSpikes = histcounts(c.poshd(cKeepIdx), hdBins);
hdRates = hdSpikes ./ (hdOccupancy .* dt);

% convert coordinates from polar to rectangular form
hdBins = -pi : 2*pi/nHdBins : pi - pi/nHdBins;
x = hdRates .* cos(hdBins);
y = hdRates .* sin(hdBins);
sumLength = sum(sqrt(x.^2 + y.^2));

% computes sum vector
rx = sum(x);
ry = sum(y);
rLength = sqrt(rx^2 + ry^2);

hdScore = rLength / sumLength;

end