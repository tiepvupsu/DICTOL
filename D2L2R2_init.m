function [D, X] = D2L2R2_init(Y, Y_range, D_range, opts)

    C = numel(Y_range) - 1;
    D = zeros(size(Y, 1), D_range(end));
    X = zeros(D_range(end), Y_range(end));
    for c = 1: C 
        range_Yc = get_range(Y_range, c);
        range_Dc = get_range(D_range, c);
        Yc = Y(:, range_Yc);
        [Dc, Xcc] = ODL(Yc, D_range(c+1) - D_range(c), opts.lambda1, opts);
        D(:, range_Dc) = Dc;
        X(range_Dc, range_Yc) = Xcc;
    end 
end 