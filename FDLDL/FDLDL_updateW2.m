function [W, iter] = FDLDL_updateW2(Y, D, X, W, Y_range, opts)
%     warning off;
    if nargin == 0 
        addpath('../utils')
        addpath('../dict_learn')
        addpath('../myQP');
        profile off;
        profile on;
        C = 100; k = 7; d = 500; N = 7;
        C = 3; k = 3; d = 3; N = 3;
        opts.K = C*k;
        Y = normc(rand(d, N*C));
        D = normc(rand(d, k*C));
        X = randn(size(D, 2), size(Y, 2));
%         W = rand(k*C, C);
        Y_range = N*(0:C);
        
        opts.lambda1 = 0.1;
        opts.lambda2 = 0.1;
        opts.lambda3 = 0.002;
        opts.delta = k;
        opts.checkgrad = 0;
        opts.showprogress = 1;
        opts.max_iter = 30;
        opts.tol = 1e-6;
        W = ones(size(D,2), C)/C;
        check_cost = true;
    end 
    if ~exist('check_cost', 'var')
        check_cost = false;
    end 
    %%
    K = opts.K;
    C = numel(Y_range) - 1;   
    %%
    X2 = zeros(size(X));
    for c = 1: C 
        range_c = Y_range(c)+1:Y_range(c+1);
        wc = W(:,c);
        X2(:, range_c) = myMult_diag_a_B(wc, X(:, range_c));
    end 
    %%
    DtD = D'*D;
    DtD2 = DtD.^2;
    H = DtD.*(X*X');
    Dhat = 2*opts.lambda3*diagonal_zero(DtD2);
    %%
    if nargin == 0
        BigD = zeros(K*C, K*C);    
        D_range = K*(0:C);
        for c = 1: C 
            for j = 1: C 
                if j == c 
                    continue 
                end 
                range_c = D_range(c) + 1: D_range(c+1);
                range_j = D_range(j) + 1: D_range(j+1);
                BigD(range_c, range_j) = Dhat;
            end 
        end 
        BigD = BigD;
        %%
        BigH = zeros(K*C, K*C);
        for c = 1: C         
            range_c = D_range(c) + 1: D_range(c+1);
            BigH(range_c, range_c) = H;        
        end 
    end
    %%
    B = zeros(K, C);    
    DtY = D'*Y;
    for c = 1: C 
        DtYc = get_block_col(DtY, c, Y_range);
        Xc = get_block_col(X, c, Y_range);
        B(:, c) = sum(DtYc .* Xc, 2);
    end 
    B = -B;
    %%
    function cost = calc_cost(W)
        w = vec(W);
        cost = 0.5*w'*(BigH + BigD)*w + vec(B)'*w;       
    end 
    %% verify cost 
    if check_cost
        W1 = rand(size(W)); cost1 = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        cost2 = calc_cost(W1);
        dif1 = cost1 - cost2;
        
        W1 = rand(size(W)); cost1 = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        cost2 = calc_cost(W1);
        dif2 = cost1 - cost2;
        
        W1 = rand(size(W)); cost1 = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        cost2 = calc_cost(W1);
        dif3 = cost1 - cost2;
        
        disp([dif1 dif2 dif3]);
    end
    %%
    if nargin == 0
        cost_slow = calc_cost(W)
    end      
    [W, iter] = solveNQP_W_ADMM_mat2(H - Dhat, B, W);    
    if nargin == 0         
        cost_slow = calc_cost(W)           
        pause
    end
end 