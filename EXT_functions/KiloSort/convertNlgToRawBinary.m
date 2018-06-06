function ops = convertNlgToRawBinary(ops)
dbstop if error;

fclose('all');

% check if file already exists
if exist(ops.fbinary, 'file')
   
    prompt = sprintf('%s already exist. Overwrite? Y/N [Y]: ', ops.fbinary);
    
    str = input(prompt, 's'); if isempty(str), str = 'Y'; end
    
    if strcmpi(str, 'y'), delete(ops.fbinary);
    elseif strcmpi(str, 'n'),  return; end
    
end

files = subdir(fullfile(ops.root, sprintf('*.dat'))); % get list of files to process
dataOut = cell(length(files), 1);
line = 0;
for ii = 1:length(files)
    fprintf(repmat('\b', 1, line));
    line = fprintf('%i/%i ,%.2f %.2f', ii, length(files), ii*100/length(files), toc);
    fidIn = memmapfile(files(ii).name,... % memory-map file for quicker processing
        'Format', 'int16');
    dataDouble = -double(fidIn.data); % convert to double and invert
    dataReshape = reshape(dataDouble, ops.NchanTOT, length(fidIn.data)/ops.NchanTOT); % reshape to 16-by-inf  
    dataOut{ii} = dataReshape'; % append to matrix (very inefficient)
end
dataOut = cell2mat(dataOut);
fidOut = fopen(ops.fbinary, 'a'); % open file for appending
fwrite(fidOut, dataOut, 'int16'); % append to binary file
fclose(fidOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% for ii = 1:ops.NchanTOT
%     line1 = fprintf('Channel %i/%i\t', ii, ops.NchanTOT);
% 
%     for kk = 1:length(files)
%         line2 = fprintf('%i/%i %.2f%%', kk, length(files), kk*100/length(files));
%     
%         channelSize = files(kk).bytes/ops.NchanTOT; % find size of each channel IN BYTES
%         offset = channelSize*(ii-1); % set offset of specific channel inside each file
%         
%         % memmap the file starting at OFFSET
%         fid = memmapfile(files(kk).name,...
%             'Format', 'int16',...
%             'Offset', offset);
%         
%         % get the block size IN SAMPLES
%         if ii == 1 && kk == 1
%             blockSize = length(fid.Data)/16;
%         end
%         
%         % add to cell array (samples x channel)
%         allData{ii, kk} = fid.Data(1 : blockSize)';   
%         
%         fprintf(repmat('\b', 1, line2));
%     end
%     fprintf(repmat('\b', 1, line1));
% end
% 
% % convert to microvolt
% Data = (-double(cell2mat(allData)) - 2048).*3.3;