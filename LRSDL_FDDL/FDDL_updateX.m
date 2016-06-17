function X = FDDL_updateX(Y, Y_range, D, D_range, X, opts)
    % X = argmin_X 0.5\|Yhat - Dhat X\| + 
    %    + 0.5*lambda2(\sum (normF2(Xi - Mi) - normF2(Mi - M)) + normF2(X) + normF2(X0 - M0))
    %    + lambda1*norm1(X)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: test module ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0     
        addpath('../utils');
        addpath('../sparse_coding');
        tic
        d       = 30;
        N       = 7;
        k       = 5;
        C       = 3 ;
        Y       = normc(rand(d,N*C));
        D       = normc(rand(d,k*C));
        Y_range = N* (0:C);
        D_range = k* (0:C);
        X       = randn(size(D,2), size(Y,2));

        opts.k0       = k0;
        opts.lambda1  = 0.01;
        opts.lambda2  = 0.002;
        opts.lambda3  = 0.1;
        opts.max_iter = 250;
        opts.show     = true;        
        opts.check_grad = true;
        opts          = initOpts(opts); % other attributes
    end
    
    lambda1  = opts.lambda1;
    lambda2  = opts.lambda2;
    DtD      = D'*D;
    D_0      = buildMhat(DtD, D_range, D_range);
    Dhat     = D_0 + 2*opts.lambda2*eye(size(D_0,1));
    DtY     = D'*Y;
    Y_0     = buildMhat(DtY, D_range, Y_range);
    %% cost w.r.t. X, X0, not includeing norm1 terms
    function cost = calc_f(X)        
        cost = 0.5*(normF2(Y - D*X) + ...
                    FDDL_fidelity(Y, Y_range, D, D_range, X)) + ...       
               0.5*opts.lambda2* (FDDL_discriminative(X, Y_range));
    end 
    %% Total cose 
    function cost = calc_F(X)
        cost = calc_f(X) + lambda1*norm1(X);
    end 
    %%
    %% Gradient for FISTA 
    function g = grad(X)
        g       = Dhat*X - Y_0 + buildM_2Mbar(X, Y_range, lambda2);
    end
    %% check gradient
    if opts.check_grad &&~check_grad(@calc_f, @grad, X)
        fprintf('Check gradient or cost again!\n')
        pause
    end       
    %% ========= Main FISTA ==============================
    optsXX0          = opts;
    optsXX0.max_iter = 300;
    L = max(eig(Dhat)) + 6*lambda2;  
    X      = fista(@grad, X, L, opts.lambda1, opts, @calc_F);
    %%
    if nargin == 0   
        fprintf('done, press any key to see results\n');
        pause;
    end
end 

