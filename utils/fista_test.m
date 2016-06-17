function FISTA_test()
    warning off all;
    addpath(fullfile('..', 'utils'));
    addpath(fullfile('..', 'build'));
    %% ========= Test mode (if nargin == 0) ==============================
    d      = 5; % data dimension
    N      = 10; % number of samples 
    k      = 10; % dictionary size 
    lambda = 0.1;
    Y      = normc(rand(d, N));
    D      = normc(rand(d, k));
    Xinit  = zeros(size(D,2), size(Y, 2));
    %
    opts.max_iter     = 300;
    opts.show_progress = 0;
    opts.checkgrad    = 0;  
    % param for mex
    param.lambda     = lambda; % not more than 20 non-zeros coefficients
    param.lambda2    = 0;
    param.numThreads = 1; % number of processors/cores to use; the default choice is -1
    param.mode       = 2;        % penalized formulation
    % mex solution and optimal value 
    addpath(fullfile('..', 'build_spams' ));
    X0       = mexLasso(Y, D, param); 
    costmex = calc_F(X0)  
    %%
    %% 
    function cost = calc_f(X)
        cost = 1/2 *normF2(Y - D*X);
    end 
    %% cost function 
    function cost = calc_F(X) 
        cost = calc_f(X) + lambda*norm1(X);
    end 
    %%
    DtD = D'*D;
    DtY = D'*Y;
    function res = grad(X) 
        res = DtD*X - DtY;
    end 
    %% Checking gradient 
    check_grad(@calc_f, @grad, Xinit);
    %%
    L = max(eig(DtD));
    opts.tol = 1e-10;
    opts.max_iter = 500;
    [X] = fista(@grad, Xinit, L, lambda, opts);
    norm(X - X0)
    costfista = calc_F(X);
    costfista - costmex
end 