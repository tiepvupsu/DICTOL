function D = DLSI_updateD(D, E, F, A, lambda, opts)
% function D = DLSI_updateD(D, E, F, A, lambda, opts)
% problem: `D = argmin_D -2trace(ED') + trace(FD'*D) + lambda *||A*D||F^2,` 
% subject to: `||d_i||_2^2 <= 1`
% where F is a positive semidefinite matrix
% ========= aproach: ADMM ==============================    
% rewrite: `[D, Z] = argmin -2trace(ED') + trace(FD'*D) + lambda ||A*Z||_F^2,` 
%     subject to `D = Z; ||d_i||_2^2 <= 1`
% aproach 1: ADMM.
% 1. D = -2trace(ED') + trace(FD'*D) + rho/2 ||D - Z + U||_F^2, 
%     s.t. ||d_i||_2^2 <= 1
% 2. Z = argmin lambda*||A*Z|| + rho/2||D - Z + U||_F^2
% 3. U = U + D - Z
% solve 1: D = argmin -2trace(ED') + trace(FD'*D) + rho/2 ||D - W||_F^2 
%                       with W = Z - U;
%            = argmin -2trace((E - rho/2*W)*D') + 
%               trace((F + rho/2 * eye())*D'D)
% solve 2: derivetaive: 0 = 2A'AZ + rho (Z - V) with V = D + U 
% `Z = B*rhoV` with `B = (2*lambda*A'*A + rho I)^{-1}`
% `U = U + D - Z` 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

    if nargin == 0
        clc;
        d = 300; 
        N = 10;
        k = 5;
        k2 = 495;
        load('tmp.mat');
        lambda = 0.01;        
        
        opts.show = 0;
        opts.max_iter = 300;   
        opts.verbose = 1;
       
    end 
    if nargin == 6
        opts.lambda = lambda;
    elseif nargin == 5
        opts = lambda;
    end                
    %%
    function cost = calcost(D)       
        cost =  -2*trace(E*D') + trace(F*D'*D) + lambda*normF2(A*D);       
    end 
    %%
    iter = 0;
    rho = 1.0;
    Z_old = D;
    U = zeros(size(D));
    I_k = eye(size(D,2));
    % B = inv(2*lambda*A'*A + rho*I_k2); However, this might be very expensive if size(A, 2) is big, which is common    
    % Instead, we can use the Sherman–Morrison formula at
    % https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula#Generalization_(Woodbury_Matrix_Identity)
    X = 2*lambda/rho*A';
    Y = A;
    B1 = X*inv(eye(size(Y, 1)) + Y*X);
    tol = 1e-5;
    optsD.max_iter = 100;
    optsD.tol = 1e-8;
    while iter < opts.max_iter 
        iter = iter + 1;
        %% ========= update D ==============================         
        W = Z_old - U;
        E2 = E + rho/2 * W;
        F2 = F + rho/2*I_k; 
        D = ODL_updateD(D, E2, F2, optsD);
        %% ========= update Z ==============================
        V = D + U;
        % Z_new = rho*B*V; slow 
        Z_new = rho*(V - B1*(Y*V)); % fast 
        e1 = normF2(D - Z_new);
        e2 = rho*normF2(Z_new - Z_old);
        if (e1 < optsD.tol && e2 < optsD.tol)
            break;
        end
        if opts.verbose
            cost = calcost(D);
            fprintf('iter = %3d | costD = %5.4f | normF2(D - Z) = %5.4f | rho(Z_new - Z_old) = %5.4f\n', iter, cost, e1, e2);
        end 
        %% ========= update U ==============================
        U = U + D - Z_new;
        Z_old = Z_new;
    end 
%     disp(t1)
%     disp(t2)
    if nargin == 0
        D = [];
    end 
end 
