function varargout = index_to_keep(st, p, s)
%% INDEX_TO_KEEP a function to clean positionial data structure based on throw_away_times and session limits
%   INPUT:  ST - structure to clean
%           p - contains throw_away_times
%           s - contains session data
%   OUTPUT: CLEAN STRUCT and keepIdx (if invoked)

sessionTime = find(st.timestamps <= s.start_time | st.timestamps >= s.end_time); % get only session times
% speed_threshold = find(st.speed >= 50 | st.speed <= 2); % get speed threshold
noBehavior = find(PRE_throw_away_times(st.timestamps, p.throw_away_times)); % remove throw away times
throwIdx = unique([sessionTime noBehavior]);

if isempty(throwIdx), keepIdx = 1:length(st.timestamps);
else, keepIdx = setdiff(1:length(st.timestamps), throwIdx);
end

varargout{1} = structfun(@(x) x(keepIdx), st, 'UniformOutput', false);

if nargout == 2
    varargout{2} = keepIdx;
end

end