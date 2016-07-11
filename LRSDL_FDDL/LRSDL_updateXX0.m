function [X, X0] = LRSDL_updateXX0(Y, Y_range, D, D_range, D0, X, X0, opts)
    if nargin == 0     
        addpath('../utils');
        addpath('../sparse_coding');
        tic
        d       = 30;
        N       = 7;
        k       = 5;
        k0      = 5;
        C       = 3 ;
        
        Y       = normc(rand(d,N*C));
        D       = normc(rand(d,k*C));
        D0      = normc(rand(d, k0));
        Y_range = N* (0:C);
        D_range = k* (0:C);
        X       = randn(size(D,2), size(Y,2));
        X0      = randn(size(D0, 2), size(Y, 2));
        
        opts.k0       = k0;
        opts.lambda1  = 0.01;
        opts.lambda2  = 0.002;
        opts.lambda3  = 0.1;
        opts.max_iter = 250;
        opts.show     = true;    
        opts.verbose   = true;
        opts          = initOpts(opts); % other attributes
    end
    %% 
    function cost = calc_f(X1)
        [X, X0] = extractFromX1(X1);
        Ybar = Y - D0*X0;
        cost = 0.5*(normF2(Ybar - D*X) + ...
                    FDDL_fidelity(Ybar, Y_range, D, D_range, X)) + ...       
               0.5*opts.lambda2* (FDDL_discriminative(X, Y_range) + ...
                                    normF2(X0 - buildMean(X0)));
    end 
    %% Total cose 
    function cost = calc_F(X1)
        cost = calc_f(X1) + lambda1*norm1(X1);
    end 
    %% [X, X0] = [X; X0]
    function [X, X0] = extractFromX1(X1)
        X = X1(1:D_range(end), :);
        X0 = X1(D_range(end) + 1: end, :);
    end    
    %% Gradient for FISTA 
    function g1 = grad(X1)
        [X, X0] = extractFromX1(X1);
        DtY     = DtY0 - DtD0*X0;
        Y_0     = buildMhat(DtY, D_range, Y_range);
        g       = Dhat*X - Y_0 + buildM_2Mbar(X, Y_range, lambda2);
        g0      = A*X0 - D0tY2 + D0tD*buildMhat(X, D_range, Y_range)...
                      - lambda2*buildMean(X0);        
        g1      = [g; g0];
    end
    %% 
    if opts.k0 == 0 
        X = FDDL_updateX(Y, Y_range, D, D_range, X, opts);
        X0 = [];
    else 
        %% Prepare data 
        lambda1  = opts.lambda1;
        lambda2  = opts.lambda2;
        DtD      = D'*D;
        D_0      = buildMhat(DtD, D_range, D_range);
        Dhat     = D_0 + 2*opts.lambda2*eye(size(D_0,1));
        D0tD0    = D0'*D0;
        A        = 2*D0tD0 + opts.lambda2*eye(size(D0,2));
        DtY0     = D'*Y;
        DtD0     = D'*D0;
        D0tD     = DtD0';
        D0tY2    = 2*D0'*Y;
        %% check gradient
        if opts.check_grad &&~check_grad(@calc_f, @grad, [X; X0])
            fprintf('Check gradient or cost again!\n')
            pause
        end       
        %% ========= Main FISTA ==============================
        optsXX0          = opts;
        optsXX0.verbose = false;
        optsXX0.max_iter = 300;
        L = max(eig(Dhat)) + max(eig(A)) + 4*lambda2 + 1;  
        X1 = [X; X0];
        X1      = fista(@grad, X1, L, opts.lambda1, optsXX0, @calc_F);
        [X, X0] = extractFromX1(X1);
    end 
    %%
    if nargin == 0   
        fprintf('done, press any key to see results\n');
        pause;
    end
end 

