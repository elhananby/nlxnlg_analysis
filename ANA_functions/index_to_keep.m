function keepIdx = index_to_keep(arr, p, s)

session_time = arr.timestamps <= s.start_time & arr.timestamps >= s.end_time; % get only session times
speed_threshold = arr.speed >= 50 | arr.speed <= 1; % get speed threshold
no_behavior = PRE_throw_away_times(arr.timestamps, p.throw_away_times); % remove throw away times
keepIdx = ~(session_time | speed_threshold | no_behavior); % NOT throw away to get indices to keep

end