function cost = FDLDL_cost(Y, D, X, W, Y_range, opts)
    if nargin == 0
        addpath('../utils');
        d = 3;
        C = 10;
        n = 30;
        k = 7;
        N = n*C;
        K = k*C;
        
        Y = normc(rand(d, N));
        D = normc(rand(d, K));
        X = randn(K, N);
        
        W = rand(K, C);
        Y_range = n*(0:C);
        opts.lambda1 = 0.001;
        opts.lambda2 = 0.01;
        opts.lambda3 = 0.001;
    end 
    %%

    cost    = 0;
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = opts.lambda3;
    w = W(:);
    % cost = 0.5*normF2(Y - D*X) + lambda1*norm1(X) + lambda2*normF2(X) + opts.lambda3*normF2(W)/2;
    cost = lambda1*norm1(X) + 0.5*lambda2*normF2(X);
    % cost0 = 0;
    m = mean(X, 2);
    C = size(W, 2);
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
        mc   = mean(Xc, 2);
        Mc   = repmat(mc, 1, size(Yc, 2));
        M    = repmat(m, 1, size(Yc, 2));
        cost = cost + lambda2*(normF2(Xc - Mc) - normF2(Mc - M));
    end 
    
    % cost = cost*2; % similar to FDDL
    
%%%%%%% Calculate the third term by definition
%     tic 
%     cost31 = 0;
%     for j = 1: C 
%         for l = 1: C 
%             if l == j 
%                 continue;
%             end 
%             for n = 1:K
%                 for m = 1: K
%                     if m == n 
%                         continue;
%                     end 
%                     cost31 = cost31 + W(m,j) * (D(:,m)'*D(:,n))^2 *W(n, l);
%                 end 
%             end 
%         end 
%     end 
%     toc
    %%
%     cost = cost1 + lambda3*cost31;    

    %% Calculate the third term in an efficient way 
    % tic 
    DtD2 = (D'*D).^2;
    % DtD2(1:K+1:K^2) = 0; %% set diagonal of DtD2 to be zeros.
    wbar = sum(W, 2);
    W2 = (wbar*wbar' - W*W');
    % W2(1:K+1:K^2) = 0;
    W2 = diagonal_zero(W2);
    tmp = W2.*(DtD2);

    cost32 = sum(sum(tmp));
    % toc 
%     cost31 - cost32
    % pause

    cost = cost + lambda3*cost32;
    
end
