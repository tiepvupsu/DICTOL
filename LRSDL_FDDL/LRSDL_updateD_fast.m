function D = LRSDL_updateD_fast(Y, Y_range, D, D_range, D0, X, X0, opts)
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
    %%    
    if opts.k0 > 0
        Ybar = Y - D0*X0;
    else 
        Ybar = Y;
    end 
    FDDL_updateD_fast(Ybar, Y_range, D, D_range, X, opts);

end