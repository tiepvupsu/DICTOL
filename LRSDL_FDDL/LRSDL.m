function [D, D0, X, X0, CoefM, coefM0, opts, rt] = LRSDL(Y, train_label, opts)
% function [D, D0, X, X0, CoefM, coefM0, opts] = LRSDL(Y, train_label, opts)
% Main LRSDL algorithm

    if nargin == 0 
        addpath('../utils')
        addpath('../sparse_coding');
        addpath('../build');
        addpath('../dict_learn');
        addpath('../ODL');        
        C = 4;    N = 10; d = 10; k = 5; k0 = 4;         
        opts.D_range = k*(0:C);
        opts.initmode = 'normal';
        opts.k0 = k0;
        opts.lambda1 = 0.0001;
        opts.lambda2 = 0.01;
        opts.lambda3 = 0.00;
        opts.max_iter = 20;
        opts.show_progress = true;
        opts.rank_flag = 1;
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
    D_range = opts.D_range;
    %% Initialization
    optsinit = opts;
    optsinit.max_iter = 30;  optsinit.verbose = 0;
    if opts.verbose
        fprintf('Initializing....');
    end 
    [D, D0, X, X0] = LRSDL_init(Y, Y_range, D_range, optsinit);  
    if opts.verbose
        fprintf('done\n');    
    end 
    %% Options for subproblems
    optsX = opts;
    optsX.verbose = 0;
    optsX.max_iter = 300;
    optsX.show_progress = 0;
    
    optsD = opts;
    optsD.max_iter = 200;
    optsD.verbose = false;

    tol_XX0 = 1e-5;
    %% Start main loop    
    iter = 0;
    maxiter = opts.max_iter;
    if opts.verbose
        cost_old = LRSDL_cost(Y, Y_range, D, D0, D_range, X, X0, opts);
        fprintf('Initial cost: %4.4f\n', cost_old);
    end 
    tic;
    while iter < opts.max_iter
        % tic
        iter = iter + 1;
        %% ========= Update X, X0 ==============================
        if opts.verbose
            fprintf('iter %3d/%3d |', iter, maxiter);
            fprintf('updating X, X0...');
        end 
        optsX.k0 = opts.k0;
        [X, X0] = LRSDL_updateXX0(Y, Y_range, D, D_range, D0, X, X0, optsX);
        %% ========= reduce shared dictionary (if needed)=================
        g0 = sum(abs(X0), 2);
        unusedid0 = find(g0 < tol_XX0);
        if numel(unusedid0) > 0
            opts.k0 = opts.k0 - numel(unusedid0);            
            D0(:, unusedid0) = [];
            X0(unusedid0, :) = [];           
        end 
        %% ========= reduce class-specific dictionaries (if needed)=======
        % g = sum(abs(X), 2);
        % unusedid = find(g < tol_XX0);
        % if numel(unusedid) > 0            
        %     D(:, unusedid) = [];
        %     X(unusedid, :) = [];
        %     opts.D_range = range_reduce(opts.D_range, unusedid);   
        %     D_range = opts.D_range;
        % end 
        %% ========= Update D ==============================       
        if opts.verbose
            fprintf('updating D...');
        end 
        optsD.k0 = opts.k0;
        D = LRSDL_updateD_fast(Y, Y_range, D, D_range, D0, X, X0, optsD);
        %% ========= Update D0 ==============================
        if opts.k0 > 0
            if opts.verbose
                fprintf('updating D0...');
            end 
%             D0 = LRSDL_updateD0(Y, D, D0, Y_range, D_range, opts, X, X0, optsD);
            D0 = LRSDL_updateD0(Y, Y_range, D, D_range, D0, X, X0, optsD);
            if opts.verbose
                if opts.show_cost 
                    cost_new = LRSDL_cost(Y, Y_range, D, D0, D_range, X, X0, opts);
%                     t = toc;
                    fprintf('K = %3d, k0 = %3d | cost_new = %5.4f', ...
                        size(D,2), size(D0,2), cost_new);
                    if abs(cost_new - cost_old) < 1e-6
                        break;
                    end
                else 
%                     t = toc;
                    fprintf('K = %3d, k0 = %3d', size(D,2), size(D0,2));
                end
            end
        end     

        %% ========= Show learn bases ==============================
        
        if opts.showD 
            if mod(sqrt(size(D, 1)), 1) ~= 0
                fprintf('Displaying bases does not apply to non square training images');
                opts.showD = false;
            elseif opts.k0 > 0 
                display_network([D, D0]);
            else 
                display_network(D);
            end 
        end         
        %% Estimated remaining time 
        t0 = toc;
        if opts.verbose
            t = t0*(opts.max_iter - iter)/iter;
            time_estimate(t);
        end
        if t0 > 20*3600 % > 20h 
            break;
        end 
    end 
    rt = toc;
    %% Refine X, X0 one more time    
    [X, X0] = LRSDL_updateXX0(Y, Y_range, D, D_range, D0, X, X0, optsX);
    opts.D_range = D_range;
    %% Mean vectors
    if opts.k0 > 0
        coefM0 = mean(X0, 2);
    else
        coefM0 = [];
    end 
    CoefM = zeros(size(X, 1), C);
    for c = 1: C
        Xc = get_block_col(X, c, Y_range);
        CoefM(:, c) = mean(Xc,2);
    end 
end 

        









