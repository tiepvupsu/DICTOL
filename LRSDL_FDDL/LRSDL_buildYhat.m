function Yhat = LRSDL_buildYhat(Y, D, X, D_range, Y_range)
    % Yhat = [Yhat1 ... Yhatj ... YhatC]
    % Yhat_c = Y_c - Dc Xcc;
    C = numel(Y_range) - 1;
    Yhat = zeros(size(Y));
    for c = 1: C 
        Yc = get_block_col(Y, c, Y_range);
        Dc = get_block_col(D, c, D_range);
        Xcc = get_block(X, c, c, D_range, Y_range);
        Yhat(:, Y_range(c) + 1: Y_range(c+1)) = Yc-Dc*Xcc;
    end
end 