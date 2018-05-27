function idx_throw = PRE_throw_away_times(t, throw)
%%PRE_THROW_AWAY_TIMES gives a (logical) vector of indices to throw away
throw = throw .* 1e6; % convert to microseconds
[m, n] = size(throw);
idx_throw = false(1, length(t))';

for i_m = 1:m % go over defined timestamps
    if length(throw(i_m,:)) == 1 % if there is only one timestamp
        if i_m == 1 % if first timestamps
            % throw away everything before
            temp_idx_throw = t <= throw(i_m);
        elseif i_m == length(m) % if last timestamps
            % throw away everything after
            temp_idx_throw = t >= throw(i_m);
        end
    elseif length(throw(i_m,:)) == 2 % if there are two timestamps
        % throw between them
        temp_idx_throw = t >= throw(i_m, 1) & t <= throw(i_m, 2);
    end
    
    if size(temp_idx_throw) == size(idx_throw)
        idx_throw = temp_idx_throw | idx_throw;
    else
        idx_throw = temp_idx_throw' | idx_throw;
    end
end
end