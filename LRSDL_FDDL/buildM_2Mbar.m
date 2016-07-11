function g2 = buildM_2Mbar(X, Y_range, lambda2)
    g2 = zeros(size(X));
    C = numel(Y_range) - 1;
    m = mean(X, 2);

    for c = 1: C 
        Xc = get_block_col(X, c, Y_range);
        mc = mean(Xc, 2);
        g2(:, Y_range(c) + 1: Y_range(c+1)) = ....
            repmat(lambda2*(m - 2*mc), 1, Y_range(c+1) - Y_range(c));
    end 
end 