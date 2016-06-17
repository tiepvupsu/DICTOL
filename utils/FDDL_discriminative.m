function cost = FDDL_discriminative(X, Y_range)

    nClasses = numel(Y_range) - 1;
    m = mean(X,2);
    cost = normF2(X);
    for c = 1: nClasses
        Xc = get_block_col(X, c, Y_range);
        Mc = buildMean(Xc);
        N = size(Xc,2);
        M = repmat(m, 1, N);
        cost = cost + normF2(Xc - Mc) - normF2(Mc - M);
    end 

end 