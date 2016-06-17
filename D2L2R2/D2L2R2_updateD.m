function D = D2L2R2_updateD(Y, Y_range, D, D_range, X, opts)

    C = numel(Y_range) - 1;
    opts.YXt = true;
    for c = 1 : C 
        Mc  = Y  - remove_block_col(D, c, D_range)*remove_block_row(X, c, D_range);
        Yc  = get_block_col(Y, c, Y_range);
        Xc  = get_block_row(X, c, D_range);
        Xcc = get_block_col(Xc, c, Y_range);
        E   = Mc*Xc' + Yc * Xcc';
        F   = 2*Xc*Xc';

        range_Dc = get_range(D_range, c);
        Dc = D(:, range_Dc);
        D(:, range_Dc) = min_rank_dict(Dc, E, F, opts.alpha, opts);
    end 
end 