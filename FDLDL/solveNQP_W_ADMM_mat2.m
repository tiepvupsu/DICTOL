function [W, iter] = solveNQP_W_ADMM_mat2(H, P, W)
    if nargin == 0
%         profile off;
%         profile on;
        addpath('../utils');
        K = 100; C = 2;
        A = rand(K, K); H = A'*A;
        P = rand(K, C);
        W = normalizeW(rand(K, C));
        opts.checkgrad = 0;
        opts.beta = 1;
    end 
    K = size(H, 1);
    C = numel(W)/K;

    rho = 100;
    Hhat = H/rho + eye(K); Hbar = inv(Hhat);

    %% =========  MAIN ADMM ==============================
    % w = 
    % w = W(:);
    % p = P(:);
    b = ones(K, 1);
    Z = zeros(size(W));
    U = Z;
    max_iter = 1000;
    iter = 0;
    alpha = .9;
    tol = 1e-5;
    while iter < max_iter
        iter = iter + 1;
        % --------------- update x -------------------------
        V = Z-U - P/rho;
        % V = reshape(v, K, C);
        % xnu = H*([z-u - c/rho; b]);        
        % w = vec(Hbar*V)  - repmat(Hbar*mean(V, 2) - b/C, C, 1);
        HV = Hbar*V;
        W = HV  - repmat(mean(HV, 2) - b/C, 1, C);
        
        % x = xnu(1: d, :);
        % norm(A*x - b)
        W = .9*W + .1*Z;
        % --------------- update z -------------------------
        Z_new = max(0, Z + U);
%         eps1 = 
        if mod(iter, 20) == 0
            if norm(W - Z_new) < tol & rho*norm(Z_new - Z) < tol 
    %         if max(vec(abs(W-Z_new))) < tol & rho*max(vec(abs(Z_new - Z))) < tol
                break;
            end 
        end 
        Z = Z_new;
        % pause
        % --------------- update u -------------------------
        U = U + W - Z;

        
    end 
    W = max(0, W);
%     W = reshape(w, K, C)

end 