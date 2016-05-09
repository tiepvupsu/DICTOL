function [D, iter] = FDLDL_updateD(Y, D, X, W, Y_range, opts)
    if nargin == 0 
        addpath('../utils')
        addpath('../dict_learn')
%         profile off;
%         profile on;
        C = 100; k = 7; d = 500; N = 7;
        aa = 50;
        C = aa; k = aa; d = aa; N = 2*aa;
        opts.K = C*k;
        Y = normc(randn(d, N*C));
        D = normc(randn(d, k*C));
        Xinit = randn(size(D, 2), size(Y, 2));
        W = rand(k*C, C);
        Y_range = N*(0:C);
        
        opts.lambda1 = 0.1;
        opts.lambda2 = 0.1;
        opts.lambda3 = 0.0001;
        opts.checkgrad = 0;
        opts.showprogress = 1;
        opts.max_iter = 300;
        opts.tol = 1e-8;
        X = randn(opts.K, size(Y,2));
        verify_grad = false;
        opts.sol = 'true'; % other options: 'paper', 'true'
    end 
    if ~exist('verify_grad', 'var')
        verify_grad = false;
    end
    %%
    if ~isfield(opts, 'sol')
        opts.sol = 'modified';
    end
    %% 
    C = numel(Y_range) - 1;
    if isfield(opts, 'delta');
        opts.delta = round(opts.K / C);
    end
    
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = opts.lambda3;

    K = opts.K;
    C = numel(Y_range) - 1;
    d = size(D, 1);

    %%
    % X2 = zeros(size(X));
    % for c = 1: C 
    %     range_c = Y_range(c)+1:Y_range(c+1);
    %     wc = W(:,c);
    %     X2(:, range_c) = myMult_diag_a_B(wc, X(:, range_c));
    % end 
    % G = X2*X2'; % \Gamma 
    % L = Y*X2'; %Lambda


    E = 0;
    XXt = X*X';
%     ww = zeros(size(XXt));
    for c = 1: C 
        Yc = get_block_col(Y, c, Y_range);
        Xc = get_block_col(X, c, Y_range);
        wc = W(:, c);
        E = E + Yc*(myMult_diag_a_B(wc, Xc))';        
    end 
%     F = XXt + XXt.*(W*W');
    F = XXt.*(W*W');

    L = E;
    G = F;


    %%
    wbar = sum(W, 2);
    W2 = (wbar*wbar' - W*W');
%     W2(1:K+1:K^2) = 0; % set diagonal = 0 
    W2 = W2 - diag(diag(W2));
    %%
    function cost = calc_f(D)
%         
%         cost = normF2(Y - D*X2);
%         DtD2 = (D'*D).^2;
%         tmp = W2.*(DtD2);
%         cost = cost + lambda3*sum(sum(tmp));
        
        cost = FDLDL_cost(Y, D, X, W, Y_range, opts);
    end 
    %% Calculate gradient w.r.t. d_n 
    function res = gradd_n(D, n)
%         res = L(:,n) - D*G(:,n) ; res1 = -res * 2;
        res = 2*(D*X2 - Y)*X2(n, :)';
%         [res1 res]
%         pause
        wbar_n = W2(:, n);
        D2 = D.*repmat(wbar_n', d, 1);            
        Q = D2*D'; % = \su
        res = res + 4*lambda3*Q*D(:,n);
    end 
    %% num_grad_dn 
    function res = num_gradd_n(D, n)
        d = size(D, 1);
        res = zeros(d, 1);
        ep = 1e-5;
        for i = 1: d 
            Dp = D;
            Dm = D;
            Dp(i,n) = Dp(i,n) + ep; fp = calc_f(Dp);
            Dm(i,n) = Dm(i,n) - ep; fm = calc_f(Dm);
            res(i)  = (fp - fm) / (2*ep);
        end            
    end 
    %% verify gradient
    if verify_grad
        n = 1;
        g1 = gradd_n(D, n);
        g2 = num_gradd_n(D, n);
        [g1 g2]
        norm(g1 - g2)
        pause
    end 
    %%
    if nargin == 0
        cost_before = FDLDL_cost(Y, D, X, W, Y_range, opts);
        fprintf('Cost before : %5f\n', cost_before);
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
            I_Q = diag(Q);     
            switch opts.sol
                case 'paper'
                    %% paper solution, this fast but higly incorrect, sometime it
                    % increases the cost function                    
                    u = (G(n,n) + 2*lambda3*I_Q).^(-1).*L(:, n) - D*G(:, n) - 2*lambda3*(Q)*D(:,n);                     
                case 'modified'
                    % my modified solution, fast, more correct
                    u = (G(n,n) + 2*lambda3*I_Q).^(-1).*L(:, n) - D*G(:, n) + G(n,n)*D(:,n) - 2*lambda3*(Q - diag(I_Q))*D(:,n); % my modified solution            
                otherwise
                    % my solution, mathematically correct but slow because it
                    % requires matrix inversion
                    u = inv(G(n,n)*eye(d) + 2*lambda3*Q) * (L(:, n) - D*G(:, n) + G(n, n)*D(:,n)); 
            end        
            % calc_f(D)
            D(:, n) = u/max(norm(u, 2), 1);

        end 
%         calc_f(D)
%         FDLDL_cost(Y, D, X, W, Y_range, opts)
        if norm(D_old - D) < 1e-6
            break;
        end 
        D_old = D; 
        
    end
%     calc_f(D)
    %%
    
    if nargin == 0
        cost_after = FDLDL_cost(Y, D, X, W, Y_range, opts);
        fprintf('Cost after  : %5f\n', cost_after);
    end 

    %%

    if nargin == 0 
%         profile viewer;
        pause
    end
end 


