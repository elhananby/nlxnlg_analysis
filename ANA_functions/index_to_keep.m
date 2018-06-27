function keepIdx = index_to_keep(arr, p, s)

sessionTime = find(arr.timestamps <= s.start_time | arr.timestamps >= s.end_time); % get only session times
% speed_threshold = find(arr.speed >= 50 | arr.speed <= 1); % get speed threshold
noBehavior = find(PRE_throw_away_times(arr.timestamps, p.throw_away_times)); % remove throw away times
throwIdx = unique([sessionTime noBehavior]);

if isempty(throwIdx), keepIdx = 1:length(arr.timestamps);
else, keepIdx = setdiff(1:length(arr.timestamps), throwIdx);
end

end