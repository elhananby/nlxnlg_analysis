function st = permute_timestamps(t)

    isi = diff(t);
    perm = isi(randperm(length(isi)));
    st = zeros(length(t), 1);
    st(1) = t(1);
    
    for ii = 2:length(t)
       st(ii) = st(ii-1) + perm(ii-1);
    end
    
end