function [X, iter] = LDL_UpdateX(Y, D, Xinit, W, Y_range, opts)
    
    if nargin == 0 
        
%         profile off;
%         profile on;
        addpath('../sparse_coding');
        addpath('../utils');
        C = 3; k = 30; d = 30; N = 30; opts.checkgrad = 0;
        C = 10; k = 5; d = 5; N = 5; opts.checkgrad = 1;
        opts.K = C*k;
        Y = normc(rand(d, N*C));
        D = normc(rand(d, k*C));
        Xinit = zeros(size(D, 2), size(Y, 2));
        % W = zeros(k*C, C);
        
        W = rand(k*C, C);      


        Y_range = N*(0:C);
        
        opts.lambda1 = 0.1;
        opts.lambda2 = 0.1;
        opts.lambda3 = 0.01;
        opts.showprogress = 0;
        opts.max_iter = 30;
        opts.tol = 1e-6;

    end 
%     lambda1 = opts.lambda1;
%     lambda2 = opts.lambda2;
%     lambda3 = opts.lambda3;
    X = Xinit;
    K = opts.K;
    C = size(W, 2);
    iter = 0;
    %%
    optsXc = opts;
    optsXc.max_iter = 300;
%     cost = LDL_cost(Y, D, X, W, Y_range, opts)
    for c = 1 : C 
        range_c = Y_range(c) + 1: Y_range(c + 1);
        Yc = Y(:, range_c);
        Xc = X(:, range_c);
        wc = W(:, c);
        Dc = myMult_A_diag_b(D, wc);

        [X(:, range_c), it] = LDL_updateXc(Yc, Dc, Xc, optsXc);
        iter = iter + it;
%         cost = LDL_cost(Y, D, X, W, Y_range, opts)
    end 
    iter = round(iter/ C);
%     pause
end 

function [X, iter] = LDL_updateXc(Y, D, Xinit, opts)
    %% solve the problem: 
    % X = \arg\min_ \|Y - DX\|_F^2 + lambda1\|X\|_1 + lambda2 \|X - M\|_F^2
    % where M = repmat(mean(X,2), 1, size(X, 2)), which is the mean matrix of X 

    lambda1 = opts.lambda1/2;
    lambda2 = opts.lambda2;
    function cost = calc_f(X)
        M = repmat(mean(X,2), 1, size(X, 2));
        cost = .5*normF2(Y - D*X) + .5*lambda2*normF2(X - M) ;
    end 
    %% cost function 
    function cost = calc_F(X) 
        cost = calc_f(X) + lambda1*norm1(X);
    end 
    %%
    DtD = D'*D;
    DtY = D'*Y;
    bigD = DtD + lambda2*eye(size(DtD,1));
    function res = grad(X) 
        M = repmat(mean(X,2), 1, size(X, 2));
        res = bigD*X - DtY - lambda2*M;
    end 
    %% Checking gradient 
%     check_grad(@calc_f, @grad, Xinit);
    %%
    L = max(eig(bigD)) + lambda2;
    opts.tol = 1e-10;
    % opts.max_iter = 300;
    [X, iter] = FISTA(@grad, Xinit, L, lambda1, opts, @calc_F);     
end 