function [W, iter] = FDLDL_updateW(Y, D, X, W, Y_range, opts)
    warning off;
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
        opts.lambda3 = 0.001;
        opts.delta = k;
        opts.checkgrad = 0;
        opts.showprogress = 1;
        opts.max_iter = 30;
        opts.tol = 1e-6;
        W = opts.delta/size(D,2) * ones(size(D,2), C);
        check_cost = true;
    end 
    if ~exist('check_cost', 'var')
        check_cost = false;
    end 
    %%
    lambda3 = opts.lambda3;
    delta   = opts.delta;
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
    function cost = calc_cost_slow(W)
        cost = 0; 
        for c = 1: C 
            Yc   = get_block_col(Y, c, Y_range);
            Xc   = get_block_col(X, c, Y_range);
            wc   = W(:, c);
            
            D_c  = myMult_A_diag_b(D, wc);
            cost = cost + .5*normF2(Yc - D_c*Xc);
            cost = cost;
            for i = 1: C 
                if i == c 
                    continue;
                else 
                    Xi = get_block_col(X, i, Y_range);
                    cost = cost + 0.5*normF2(D_c*Xi);
                end 
            end 
        end 
        wbar = sum(W, 2);
        W2 = (wbar*wbar' - W*W');
        W2(1:K+1:K^2) = 0; % set diagonal = 0 
        tmp = W2.*(DtD2);
        cost32 = sum(sum(tmp));
        cost = cost + lambda3*cost32;
    end
    %%
    DDD = DtD2;
    DDD(1:K+1: K^2) = 0;
    DtY = D'*Y;
    function [Q, p, cost] = calc_cost_fast(W, c)
        Xc = get_block_col(X, c, Y_range);
        DtYc = get_block_col(DtY, c, Y_range);
        % Q = DtD.*(Xc*Xc');
        Q = DtD.* (X*X');
        %%
        WW = W; WW(:,c) = [];
        ww = sum(WW, 2);
        wc = W(:, c);
        b = DDD'*ww;
        
        fc = -sum(DtYc .* Xc, 2);


        p = fc + 2*lambda3*b;
        if nargout == 3 
            cost = 0.5*wc'*Q*wc + p'*wc;
            % cost = 2*cost;
        end 
    end 
    %% Verify cost wc 
    if check_cost
        W1 = rand(size(W));
        k = 1;
        cost_slow = calc_cost_slow(W1); 
        cost_slow = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        [~,~,cost_fast] = calc_cost_fast(W1, k); 
        dif1 = cost_slow - cost_fast;
        
        W1(:, k) = rand(K, 1); cost_slow = calc_cost_slow(W1); 
        cost_slow = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        [~,~,cost_fast] = calc_cost_fast(W1, k); 
        dif2 = cost_slow - cost_fast;
        
        W1(:, k) = rand(K, 1); cost_slow = calc_cost_slow(W1); 
        cost_slow = FDLDL_cost(Y, D, X, W1, Y_range, opts);
        [~,~,cost_fast] = calc_cost_fast(W1, k); 
        dif3 = cost_slow - cost_fast;
        fprintf('=============== Checking cost ======================\n');
        fprintf('The following three number should be the same:\n');
        disp([dif1 dif2 dif3]); % these three value should be the same, showing that two different 
        % ways of calculating cost are similar (differences are constants)
        fprintf('====================================================\n');
        pause
    end
    %% Main algorithm
    iter = 0;
    options = optimset('Algorithm','interior-point-convex','Display','off');
    A = ones(1, K);
    if nargin == 0
        cost_slow = FDLDL_cost(Y, D, X, W, Y_range, opts)
    end
%     b = delta;
    while iter < opts.max_iter 
        iter = iter + 1;
        W_old = W;
        for c = 1: C 
           %% Each wc will be found be solving the following problem
           % wc = argmin_{wc} 0.5*wc'*W*wc + p'*wc 
           % s.t. wc >= 0 and sum(wc) = delta.
           [Q, p, cost] = calc_cost_fast(W, c);
%            cost
   
%            W(:, c) = quadprog(Q, p, [], [], ones(1, K), delta, zeros(K, 1), [], W(:, c), options);
          
           % [W(:, c), ~] = NNQP_ADMM(Q, p, A, delta, W(:,c), opts)   ;  
           W(:,c) = ConNQP(Q, p, A, -delta, W(:,c));      
           if nargin == 0
               FDLDL_cost(Y, D, X, W, Y_range, opts);
           end
        end 
%         calc_cost_slow(W)
%         cost = FDLDL_cost(Y, D, X, W, Y_range, opts)
        diff = norm(norm(W_old - W));
        if  diff < opts.tol 
            break;
        end         
    end 
%     imagesc(W);
    if nargin == 0         
        cost_slow = FDLDL_cost(Y, D, X, W, Y_range, opts)   
        
        pause
    end
end 