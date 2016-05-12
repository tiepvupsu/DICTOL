function cost = FDDL_fidelity(Y, Y_range, D, D_range, X)

    nClasses = numel(Y_range) - 1;
    cost = 0;
    for c = 1: nClasses
        Yc = get_block_col(Y, c, Y_range);
        Dc = get_block_col(D, c, D_range);
        Xc = get_block_row(X, c, D_range);
        Xcc = get_block_col(Xc, c, Y_range);
        cost = cost + normF2(Yc - Dc *Xcc);
        for j = 1:nClasses
            if j == c 
                continue;
            else
                Xcj = get_block_col(Xc, j, Y_range);
                cost = cost + normF2(Dc*Xcj);
            end 
        end 
    end 

end 