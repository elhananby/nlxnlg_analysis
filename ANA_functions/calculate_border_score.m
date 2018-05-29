


function borderScore = calculate_border_score(rateMap)
maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);

% Calculates CM
nFields = max(modifiedRateMap(:));
maxDistance = 0;
longestFieldId = 0;
maxIdx = 0;
for i = 1:nFields
    % Determines number of adjacent pixels to a wall in a firing field
    [fieldRow, fieldCol] = find(modifiedRateMap == i);
    
    nFieldPixels = length(fieldRow);
    % Fields defined as needing to have more than 6 pixels
    if (nFieldPixels <= 6)
        %         modifiedRateMap(modifiedRateMap == i) = 0;
        continue;
    end
    
    % Calculates length of pixels along each wall
    rowLeftPixels = fieldRow(fieldCol == 1); % pixels along vertical left wall
    leftDistance = abs(max(rowLeftPixels) - min(rowLeftPixels));
    if isempty(leftDistance)
        leftDistance = 0;
    end
    
    rowRightPixels = fieldRow(fieldCol == length(rateMap)); % pixels along vertical right wall
    rightDistance = abs(max(rowRightPixels) - min(rowRightPixels));
    if isempty(rightDistance)
        rightDistance = 0;
    end
    
    colTopPixels = fieldCol(fieldRow == 1); % pixels along horizontal top wall
    topDistance = abs(max(colTopPixels) - min(colTopPixels));
    if isempty(topDistance)
        topDistance = 0;
    end
    
    colBottomPixels = fieldCol(fieldRow == length(rateMap)); % pixels along horizontal bottom wall
    bottomDistance = abs(max(colBottomPixels) - min(colBottomPixels));
    if isempty(bottomDistance)
        bottomDistance = 0;
    end
    
    % Determines maximal length of wall touching the firing field
    [fieldMaxDistance, idx] = max([leftDistance rightDistance topDistance bottomDistance]);
    
    if fieldMaxDistance > maxDistance
        longestFieldId = i;
        maxIdx = idx;
        maxDistance = fieldMaxDistance;
    end
end

CM = maxDistance;

% Calculates DM
% Commented out lines used for Method 2
% [row, col] = find(modifiedRateMap == longestField);
[row, col] = find(modifiedRateMap ~= 0);
nPixels = length(row);
totalDistance = 0;
for i = 1:nPixels
    %     if maxIdx == 1 || maxIdx == 3
    %         wall = 1;
    %     else
    %         wall = length(rateMap);
    %     end
    %
    %     if maxIdx == 1 || maxIdx == 2
    %         distance = abs(col(i) - wall);
    %     else
    %         distance = abs(row(i) - wall);
    %     end
    
    distance = min([abs(col(i) - 1) abs(col(i) - length(rateMap)) ...
        abs(row(i) - 1) abs(row(i) - length(rateMap))]);
    
    totalDistance = totalDistance + distance;
end

DM = totalDistance / nPixels;
if (isnan(DM))
    DM = 0;
end

borderScore = (CM - DM)/(CM + DM);

end