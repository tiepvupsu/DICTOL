function [D, X, W, opts] = FDLDL(Y, train_label, opts)
    % Cost function 
    %     profile on;
    if nargin == 0 
        addpath('../utils')
        addpath('../sparse_coding');
        addpath('../build');


        C = 3;    N = 7;    d = 10; k = 7;
        K = C*k; % total number of bases
        
        opts.K = K; 
        opts.lambda1 = 0.0001;
        opts.lambda2 = 0.001;
        opts.lambda3 = 0.00;
        opts.delta = k+1;
        opts.max_iter = 10;
        opts.showprogress = false;
        Y = normc(rand(d, C*N));
        train_label = [];
        opts.tol = 1e-10;
        for c = 1: C
            train_label = [train_label c*ones(1, N)];
        end 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: test module ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    d = size(Y, 1);
    Y_range = LabelToRange(train_label);
    C = numel(Y_range) - 1;
    k = opts.K/C;
    if ~isfield(opts, 'delta')
        opts.delta = opts.K / C;
    end
    
    
    C = numel(Y_range) -1 ;
    D_range = k*(0:C);
    opts.delta = k;
    delta = opts.delta 
    %% init W 
    W = ones(opts.K, C)*delta/opts.K;
    % for c = 1: C 
    %     range_c = D_range(c)+1: D_range(c+1);
    %     W(range_c, c) = k/opts.delta; 
    % end 

    
    %% init D 
    D = PickDfromY(Y, Y_range, k);
    %% init X 
    X = zeros(size(D, 2), size(Y, 2));
    tic;

    % optsinit = opts;
    % optsinit.max_iter = 5;
    % tic;
    % [D, X] = myODL(Y, opts.K, opts.lambda1, optsinit, 'fista');
    % % [W, tt] = FDLDL_UpdateW(Y, D, X, W, Y_range, opts);
    % [W, tt] = FDLDL_UpdateW2(Y, D, X, W, Y_range, opts);
    % [D, X, W, K] = reduce(D, X, W, 'W');
   
        %%
    costinit = FDLDL_cost(Y, D, X, W, Y_range, opts);
    t = toc;
    fprintf('costinit = %5f, running time: %d\n', costinit, t);
    iter = 0;
    %% 
    optsX = opts;
    optsX.max_iter = 300;
    optsD = opts;
    optsD.max_iter = 100;
    optsD.sol = 'modified';
    % optsD.sol = 'paper';
    % tic ;
    w_threshold = 1e-4;
    while iter < opts.max_iter
        iter = iter + 1;
%         disp(opts.K);
        %% ========= Update X ==============================
        [X, tt] = FDLDL_updateX(Y, D, X, W, Y_range, optsX);
        costX = FDLDL_cost(Y, D, X, W, Y_range, opts);
        t = toc;
        fprintf('iter = %3d, costX = %5f, running time: %ds, iters: %4d\n', iter, costX, t, tt);
        % reduce 
        [D, X, W, K] = reduce(D, X, W, 'X', w_threshold);
        opts.K = K;
        optsX.K = K;
        %% ========= Update D ==============================
        optsD.K = opts.K;
        [D, tt] = FDLDL_updateD(Y, D, X, W, Y_range, optsD);
        costD = FDLDL_cost(Y, D, X, W, Y_range, opts);
        if costD >= costX - 0.1 && ~strcmp(optsD.sol, 'true')
            optsD.sol = 'true'
            [D, tt] = FDLDL_updateD(Y, D, X, W, Y_range, optsD);
            costD = FDLDL_cost(Y, D, X, W, Y_range, opts);
        end
        t = toc;
        fprintf('iter = %3d, costD = %5f, running time: %ds, iters: %4d\n', iter, costD, t, tt);
        %% ========= Update W if costD is close to costX==============================
        if abs(costD - costX) < 1e-1 
            [W, tt] = FDLDL_updateW2(Y, D, X, W, Y_range, opts);
            costW = FDLDL_cost(Y, D, X, W, Y_range, opts);
            t = toc;
            fprintf('iter = %3d, costW = %5f, running time: %ds, iters: %4d\n', iter, costW, t, tt);
            [D, X, W, K] = reduce(D, X, W, 'W', w_threshold);
            opts.K = K;
            optsD.K = K;
            optsX.K = K;
            if abs(costW - costD) < 1e-5
                break;
            end 
        end 
        % imagesc(W); drawnow();       
        
    end 

    if nargin == 0
        W
        pause;
    end
end 

function [D, X, W, K] = reduce(D, X, W, which_var, w_threshold)
%         %% ========= Reduce ==============================
    if strcmp(which_var, 'X') == 1 
        X1 = sum(abs(X), 2);
    else 
        X1 = sum(abs(W), 2);
    end
    null_id = find(X1 < w_threshold);
    X(null_id,:) = [];
    D(:, null_id) = [];
    W(null_id, :) = [];
    K = numel(X1) - numel(null_id);
   
end 

