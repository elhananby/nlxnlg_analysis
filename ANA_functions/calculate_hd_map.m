function [hdOccupancy, hdSpikes, hdRates, hdScore] = calculate_hd_map(c, vt)
global nHdBins dt;

hdBins = -pi : 2*pi/nHdBins : pi;

hdOccupancy = histcounts(vt.poshd, hdBins);
hdSpikes = histcounts(c.poshd, hdBins);
hdRates = hdSpikes ./ (hdOccupancy .* dt);

hdBinsCalc = linspace(-pi, pi, nHdBins);
hdScore = circ_r(hdBinsCalc', hdRates');

% % convert coordinates from polar to rectangular form
% hdBins = -pi : 2*pi/nHdBins : pi - pi/nHdBins;
% x = hdRates .* cos(hdBins);
% y = hdRates .* sin(hdBins);
% sumLength = sum(sqrt(x.^2 + y.^2));
% 
% % computes sum vector
% rx = sum(x);
% ry = sum(y);
% rLength = sqrt(rx^2 + ry^2);

% hdScore = rLength / sumLength;
end