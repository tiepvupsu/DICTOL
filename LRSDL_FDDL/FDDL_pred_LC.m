function pred = FDDL_pred_LC(Y, D, CoefM, opts) % LC
%     vgamma = opts.gamma;
    gamma1 = opts.gamma1;
    gamma2 = opts.gamma2;
    C = size(CoefM, 2);
    E = zeros(C, size(Y, 2));
    for c = 1: C 
        Dc = get_block_col(D, c, opts.D_range);
        mii = get_block_row(CoefM(:, c), c, opts.D_range);
        [~, min_cost] = sparse_coding(Y, Dc, [], mii, gamma1, gamma2);
        E(c, :)  = min_cost;
    end 
    [~, pred] = min(E);
    % opts.max_iter = 100;
    % [X, ~] = lasso_fista(Y, D, zeros(size(D,2), size(Y,2)), vgamma, opts);
    % C = size(CoefM,2);
    % % w = 0.5;
    % E = zeros(C, size(Y,2));
    % for c = 1: C 
    %     Dc = get_block_col(D, c, opts.D_range);
    %     Xc = get_block_row(X, c, opts.D_range);
    %     R1 = Y - Dc*Xc;
    %     E1 = sum(R1.^2);
    %     R2 = X - repmat(CoefM(:, c), 1, size(Y,2));
    %     E2 = sum(R2.^2);
    %     E(c,:) = E1 + opts.weight*E2;
    % end 
    % [~, pred] = min(E);

    % sparse_coding();
end 

function [X, min_cost] = sparse_coding(Y, D, Xinit, m, gamma1, gamma2)
% function X = sparse_coding(Y, D, Xinit, m, gamma1, gamma2)
% Solve the problem:
% X = argmin_X .5*||Y - DX|| + gamma1*||X|| + 0.5*gamma2* ||X - M||_F^2 
% gradient: D'(DX - Y) + gamm2*(X - M) = (D'D + gamma2*I)*X - D'Y - gamma2*M;
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/23/2016 11:00:18 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0 %% test mode 
        d = 10;
        N = 20;
        k = 10;
        Y = normc(rand(d, N));
        D = normc(rand(d, k));
        Xinit = [];
        m = rand(k, 1);
        gamma1 = 0.001;
        gamma2 = 0.01;
    end 
    %%
    if numel(Xinit) == 0
        Xinit = zeros(size(D, 2), size(Y, 2));
    end 
    M = repmat(m, 1, size(Y, 2));
    %% cost 
    function cost = calc_f(X)
        cost = 0.5*normF2(Y - D*X) + 0.5*gamma2*normF2(X - M);
    end 
    %% cost overall 
    function cost = calc_F(X) 
        cost = calc_f(X) + gamma1*norm1(X);
    end 
    %% gradient 
    A = D'*D + gamma2*eye(size(D, 2));
    B = D'*Y + gamma2*M;
    function g = grad(X)
        g = A*X - B;
    end 
    %% check_grad 
    %  check_grad(@calc_f, @grad, Xinit);
    %% L 
    L = max(eig(A));
    %% opts 
    opts.max_iter = 300;
    opts.verbose = 0;
    [X, ~] = fista(@grad, Xinit, L, gamma1, opts, @calc_F);
    %% min_cost for each sample (column)
    min_cost = zeros(1, size(Y, 2));
    R1 = Y - D*X;
    R2 = X - M;
    min_cost = sum(R1.^2, 1) + gamma1*sum(abs(X), 1) + ...
                gamma2*sum(R2.^2, 1);
end 
