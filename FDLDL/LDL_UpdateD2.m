function [D, iter] = LDL_UpdateD2(Y, D, X, W, Y_range, opts)
    if nargin == 0 
        addpath('../utils')
        addpath('../dict_learn')
        profile off;
        profile on;
        C = 100; k = 7; d = 500; N = 7;
        aa = 10;
        C = aa; k = aa; d = aa; N = aa;
        opts.K = C*k;
        Y = normc(randn(d, N*C));
        D = normc(randn(d, k*C));
        Xinit = randn(size(D, 2), size(Y, 2));
        W = rand(k*C, C);
        Y_range = N*(0:C);
        
        opts.lambda1 = 0.1;
        opts.lambda2 = 0.1;
        opts.lambda3 = 1.00;
        opts.checkgrad = 0;
        opts.showprogress = 1;
        opts.max_iter = 30;
        opts.tol = 1e-6;
        X = rand(opts.K, size(Y,2));
        
    end 
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = opts.lambda3;

    K = opts.K;
    C = numel(Y_range) - 1;
    d = size(D, 1);

    %%
    X2 = zeros(size(X));
    for c = 1: C 
        range_c = Y_range(c)+1:Y_range(c+1);
        wc = W(:,c);
        X2(:, range_c) = myMult_diag_a_B(wc, X(:, range_c));
    end 
    G = X2*X2'; % \Gamma 
    L = Y*X2'; %Lambda
    %%
    wbar = sum(W, 2);
    W2 = (wbar*wbar' - W*W');
%     W2(1:K+1:K^2) = 0; % set diagonal = 0 
    W2 = W2 - diag(diag(W2));
    %%
    function cost = calc_f(D)
        cost = normF2(Y - D*X2);
        DtD2 = (D'*D).^2;
        tmp = W2.*(DtD2);
        cost = cost + lambda3*sum(sum(tmp));
    end 
    %%
    iter = 0;
    D_old = D;
%     calc_f(D)
    while (iter < opts.max_iter)
        iter = iter + 1;
        for n = 1: opts.K 
            %% see equation (12), (13) in 
            % http://www.cv-foundation.org/openaccess/content_cvpr_2014/papers/Yang_Latent_Dictionary_Learning_2014_CVPR_paper.pdf 
            wbar_n = W2(:, n);
%             for tt = 1: 100 
            D2 = D.*repmat(wbar_n', d, 1);            
            Q = D2*D'; % = \sum_{m\neq n} d_m d_m^T \sum_j=1}^C\sum_{l\neqj} w_{j,m}w_{l,n}
            %% paper solution
            I_Q = diag(Q);
            z = (G(n,n) - 2*lambda3*I_Q).^(-1);
            % v = L(:, n) - D*G(:,n) - 2*lambda3*Q*D(:, n);
%             v = L(:, n) - D*G(:,n) + D(:, n)*G(n,n);
            
            v = L(:, n) - D*G(:, n) + G(n,n)*D(:,n) + 2*lambda3*(Q - diag(I_Q))*D(:,n);
            u = z.*v;
            %% my solution
            u = inv(G(n,n)*eye(d) - 2*lambda3*Q) * (L(:, n) - D*G(:, n) + D(:, n)*G(n,n));
%             calc_f(D)
            D(:, n) = u/max(norm(u, 2), 1);
%             end 
%             calc_f(D)
        end 
        calc_f(D)
%         LDL_cost(Y, D, X, W, Y_range, opts)
        if norm(D_old - D) < 1e-6
            break;
        end 
        D_old = D; 
        
    end
%     calc_f(D)
    

    %%

    if nargin == 0 
        profile viewer;
        pause
    end
end 


