%  create a channel map file

Nchannels = 16;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]; % grouping of channels (i.e. tetrode groups)

fs = 312500; % sampling frequency
save('D:\new_exp_data\nlg_animal20_Day30618_1\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')