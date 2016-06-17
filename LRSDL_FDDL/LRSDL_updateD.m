function D = LRSDL_updateD(Y, D, D0, Y_range, D_range, opts, X, X0, pars)
    if nargin == 0     
        addpath('../utils');
        addpath('../sparse_coding');
        tic
        d = 30;
        N = 7;
        k = 5;
        k0 = 5;
        C = 3 ;
        Y = normc(rand(d,N*C));
        D = normc(rand(d,k*C));
        D0 = normc(rand(d, k0));
        Y_range = N* (0:C);
        D_range = k* (0:C);
        opts.lambda1 = 0.01;
        opts.lambda2 = 0.002;
        opts.lambda3 = 0.1;
        opts.max_iter = 250;
        opts.show = true;
        X = randn(size(D,2), size(Y,2));
        X0 = randn(size(D0, 2), size(Y, 2));
    end

    Z = Y - D0*X0;
    W = Z - D*X;

    %% ========= update D1 ==============================

    X1 = X(1:D_range(2), :);
    X11 = X1(:, 1: Y_range(2));
    D1 = D(:, 1: D_range(2));
    W = W + D1*X1;
    E = W*X1' + get_block_col(Z, 1, Y_range)*X11';
    F = 2*X1*X1';
    D_prev = updateD_EF(D1, E, F, pars.max_iter);
    D(:, 1: D_range(2)) = D_prev;

    %% ========= update Di, i > 1 ==============================
    for i = 2: numel(Y_range) - 1
        D_cur = get_block_col(D, i, D_range);
        D_com = [-D_prev, D_cur];
        X_com = X(D_range(i-1) + 1: D_range(i+1), :);
        W = W + D_com*X_com; %W = W - D_{i-1}X^{i-1} + D_iX^i
        Xi = X(D_range(i) + 1: D_range(i+1), :);
        Xii = Xi(:, Y_range(i) + 1: Y_range(i+1));
        E = W*Xi' + get_block_col(Z, i, Y_range)*Xii';
        F = 2*Xi*Xi';
        D_prev = updateD_EF(D_cur, E, F, pars.max_iter);
        D(:, D_range(i) + 1: D_range(i+1)) = D_prev;
    end
end