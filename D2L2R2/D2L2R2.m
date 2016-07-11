function [D, D_range, X, CoefM, opts, rt] = D2L2R2(Y, train_label, opts)
% function [D, D0, X, X0, CoefM, coefM0, opts] = D2L2R2(Y, train_label, opts)
% Main D2L2R2 algorithm    
    if nargin == 0 
        addpath('../utils')
        addpath('../sparse_coding');
        addpath('../build');
        addpath('../dict_learn');
        addpath('../ODL');        
        C = 4;    N = 10; d = 10; k = 5; k0 = 4;         
        opts.D_range = k*(0:C);
        opts.initmode = 'norma';
        opts.lambda1 = 0.0001;
        opts.lambda2 = 0.01;
        opts.alpha = 0.00;
        opts.max_iter = 20;
        opts.show_progress = true;
        opts.rank_flag = 1;
        opts.showD = false;
        opts = initOpts(opts);
        
        Y = normc(rand(d, C*N));
        train_label = [];
        opts.tol = 1e-10;
        for c = 1: C
            train_label = [train_label c*ones(1, N)];
        end 
     end
    %% Parameter preparation
    C = max(train_label);    
    Y_range = label_to_range(train_label);
    D_range = opts.k * (0:C);
    %% Initialization
    optsinit = opts;
    optsinit.max_iter = 30;  optsinit.verbose = 0;
    fprintf('Initializing....');
    [D, X] = D2L2R2_init(Y, Y_range, D_range, optsinit );  
    fprintf('done\n');    
    %% Options for subproblems
    optsX = opts;
    optsX.max_iter = 300;
    optsX.show_progress = 0;
    optsX.verbose = false;
    
    optsD = opts;
    optsD.max_iter = 200;
    optsD.show_cost = false;
    optsD.verbose = false;
    %% Start main loop    
    iter = 0;
    cost_old = D2L2R2_cost(Y, Y_range, D, D_range, X, opts);
    fprintf('Initial cost: %4.4f\n', cost_old);
    tic;
    while iter < opts.max_iter
        % tic
        iter = iter + 1;        
        %% ========= Update X ==============================        
        X = FDDL_updateX(Y, Y_range, D, D_range, X, optsX);
        t = toc;
        if t > 20*3600
            break;
        end 
        if opts.verbose        
            costX = D2L2R2_cost(Y, Y_range, D, D_range, X, opts); 
            fprintf('iter %3d/%d | costX = %5.5f\n', ...
                iter, opts.max_iter, costX);
        end 
        %% ========= Update D ==============================               
        D = D2L2R2_updateD(Y, Y_range, D, D_range, X, optsD);
        if opts.showD             
            display_network(D);
        end         
        if opts.verbose
            costD = D2L2R2_cost(Y, Y_range, D, D_range, X, opts);         
            %% Estimated remaining time 
            fprintf('               costD = %5.5f', costD);
            t = t*(opts.max_iter - iter)/iter;
            time_estimate(t);
        end 
        t = toc;
        if t > 20*3600
            break;
        end
    end 
    rt = toc;
    %% mean vectors
    CoefM = zeros(size(X, 1), C);
    for c = 1: C
        Xc = get_block_col(X, c, Y_range);
        CoefM(:, c) = mean(Xc,2);
    end 
end 

        









