function D = DLSI_updateD(Y, X, D, A, lambda, opts)
% function D = DLSI_updateD(Y, X, D, A, lambda, opts)
%problem: D = \arg\min_D \|Y - D*X\|_F^2 + \lambda *\|A*D\|F^2, 
% subject to: \|d_i\|_2^2 \leq 1
% ========= aproach: ADMM ==============================    
% rewrite: [D, Z] = arg\min \|Y - D*X\|_F^2 + \lambda \|A*Z\|_F^2, 
%     subject to D = Z; \|d_i\|_2^2 \leq 1
% aproach 1: ADMM.
% 1. D = \arg\min\|Y - D*X\| + \rho/2 \|D - Z + U\|_F^2, 
%     s.t. \|d_i\|_2^2 \leq 1
% 2. Z = \arg\min \lambda*\|A*Z\| + \rho/2\|D - Z + U\|_F^2
% 3. U = U + D - Z
% solve 1: D = \arg\min \|Y - D*X\|_F^2 + \rho/2 \|D - W\|_F^2 
%                       with W = Z - U;
%            = \arg\min -2trace((YX' - \rho/2*W)*D') + 
%               trace((X*X' + \rho/2 * eye())*D'D)
% solve 2: derivetaive: 0 = 2A'AZ + \rho (Z - V) with V = D + U 
% Z = B*\rho V with B = (2\lambdaA'*A + \rho I)^{-1}
% U = U + D - Z 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
        d = 30; 
        N = 10;
        k = 10;
        k2 = 20;
        if exist('dlsi_test.mat', 'file')
            load('dlsi_test.mat');
        else 
            Y = normc(rand(d, N));
            D = normc(rand(d, k));
            X = 0.01 * rand(k, N);
            A = normc(rand(k2, d));     
            save('dlsi_test.mat', 'Y', 'D', 'X', 'A');
        end 
        % load 
        lambda = 0.01;        
        opts.show = true;
        opts.max_iter = 100;
        % opts.DLSI_YXt = true; Y = Y*X'; X = X*X';
        opts.DLSI_YXt = false;
    end 
    if nargin == 6
        opts.lambda = lambda;
    elseif nargin == 5
        opts = lambda;
    end            
    lambda = opts.lambda;
    function cost = calcost(D)
        if isfield(opts, 'DLSI_YXt') == 1 && opts.DLSI_YXt
            cost = lambda*normF2(A*D) -2*trace(YXt*D') + trace(XXt*D'*D);
        else 
            cost = normF2(Y - D*X) + lambda*normF2(A*D);
        end                 
    end 

    iter = 0;
    rho = 1;
    Z_old = D;
    U = zeros(size(D));
    if (isfield(opts, 'DLSI_YXt')) && (opts.DLSI_YXt)
        YXt = Y;
        XXt = X;
    else
        YXt = Y*X';
        XXt = X*X';
    end    
    k2 = size(A,2);
    I_k2 = eye(k2);
    I_k = eye(size(D,2));
    B = inv(2*lambda*A'*A + rho*I_k2);
    tol = 1e-5;
    optsD.max_iter = 100;
    optsD.tol = 1e-8;
    while iter < opts.max_iter 
        iter = iter + 1;

        %% ========= update D ==============================
        W = Z_old - U;
        E = YXt + rho/2 * W;
        F = XXt + rho/2*I_k; 
        D = ODL_updateD(D, E, F, optsD);

        %% ========= update Z ==============================
        V = D + U;
        Z_new = B*V;

        e1 = normF2(D - Z_new);
        e2 = rho*normF2(Z_new - Z_old);
        if (e1 < tol && e2 < tol)
            break;
        end
        if opts.show
            cost = calcost(D);
            fprintf('iter = %3d || costD = %5f || normF2(D - Z) = %5f || rho(Z_new - Z_old) = %5f\n', iter, cost, e1, e2);
        end 
        % if opts.show
        %     cost = calcost(D);
        %     fprintf('iter = %3d || costD = %5f\n', iter, cost);
        % end 

        %% ========= update U ==============================
        U = U + D - Z_new;
        Z_old = Z_new;
    end 
end 
