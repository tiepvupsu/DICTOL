function [X, cost] = myLassoWeighted_fista(Y, D, W, Xinit, lambda, opts)
    %% ================== File info ==========================
    % Author            : Tiep Vu (http://www.personal.psu.edu/thv102/)
    % Time created      : 3/9/2016 6:10:31 PM
    % Last modified     : 3/9/2016 6:10:33 PM
    % Description       : myLassoWeighted_fista function solves the Lasso problem:    %       
    %           X = argmin_X 0.5*normF2(Y - D*X) + lambda*(W.*X)  (1) 
    %       with .* is the elemenwise product operator,
    %       using FISTA method in 
    %       1. A. Beck and M. Teboulle, “A fast iterative shrinkage-
    %       thresholding algorithm for linear inverse problems,�? SIAM
    %       Journal on Imaging Sciences, vol. 2, no. 1, pp. 183–202, 2009.
    %     
    %     INPUT:
    %       Y, D, W as in (1)
    %       Xinit: initial X, if not given (size(X) = [0, 0]), then X is the zero matrix with 
    %           appropriate size
    %       opts: a struct of parameters/options:
    %           - max_iter: maximum iterations, if not given, max_iter = 300. 
    %           - tol: tolerence - another stopping criteria
    %           - show_progress: show progress
    %     OUTPUT: 
    %       X: solution 
    %       cost: if nargout == 2, then cost is optimal value of (1).
    %% ================== end File info ==========================
    warning off all;
    addpath(fullfile('..', 'utils'));
    %% ========= Test mode (if nargin == 0) ==============================
    if nargin == 0
        d      = 500; % data dimension
        N      = 100; % number of samples 
        k      = 100; % dictionary size 
        lambda = 0.01;
        Y      = normc(rand(d, N));
        D      = normc(rand(d, k));
        W      = ones(size(D,2), size(Y,2)) + 10*eye(size(D,2));
        Xinit  = zeros(size(D,2), size(Y, 2));
        %
        opts.max_iter     = 300;
        opts.show_progress = 1;
        opts.checkgrad    = 0;  
        % param for mex
        param.lambda     = lambda; % not more than 20 non-zeros coefficients
        param.lambda2    = 0;
        param.numThreads = 1; % number of processors/cores to use; the default choice is -1
        param.mode       = 2;        % penalized formulation
        % mex solution and optimal value 
        addpath(fullfile('..', 'build_spams' ));
        X       = mexLassoWeighted(Y,D,W,param);
        costmex = calc_F(X)          
    end 

    %% ========= Main code ==============================
    if(sum(size(Xinit))) == 0
        Xinit = zeros(size(D, 2), size(Y,2));
    end
    opts = initOpts(opts);
    if ~isfield(opts, 'tol')
        opts.tol = 1e-7;
    end
    %%
    DtD = D'*D;
    DtY = D'*Y;
    %%
    function cost = calc_f(X)
        cost = 0.5*normF2(Y - D*X);
    end 
    %%
    function cost = calc_F(X)
        cost = calc_f(X) + lambda*norm1(W.*X);
    end 
    %%
    function g = grad(X)
        g = DtD*X - DtY;
    end 
    %%
    function g = num_grad(X)
        ep = 1e-4;
        g = zeros(size(X));
        for i = 1: size(X, 1)
            for j = 1: size(X,2)
                Xp = X;
                Xm = X;
                Xp(i,j) = Xp(i,j) + ep; fp = calc_f(Xp);
                Xm(i,j) = Xm(i,j) - ep; fm = calc_f(Xm);
                g(i,j)  = (fp - fm) / (2*ep);
            end 
        end 
    end 
    %%
    if opts.checkgrad %check grad 
        fprintf('checking gradient...\n');
        X = rand(size(X));
        diff_grad = (grad(X) - num_grad(X));        
        fprintf('...done, the difference should be small, diff = %5f\n', normF2(diff_grad)); 
        fprintf('Press any key to continue...\n');
        pause
    end 
    %%
    %% ========= Main FISTA ==============================
    L = max(eig(DtD));
    Linv = 1/L;    
    lambdaLiv = lambda*W*Linv;
    x_old = Xinit;
    y_old = Xinit;
    t_old = 1;
    iter = 0;
    if opts.show_progress
        fprintf('(max_iter = %d): 0.', opts.max_iter);
    end 
    if nargin == 0 
        cost_old =  calc_F(x_old);
    end 
    while  iter < opts.max_iter
        iter = iter + 1;
        x_new = shrinkage(y_old - Linv*grad(y_old), lambdaLiv);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
        % check stop criteria
        e = norm1(x_new - x_old)/numel(x_new);
        if e < opts.tol
            break;
        end
        % update
        x_old = x_new;
        t_old = t_new;
        y_old = y_new;
        % show progress
        if opts.show_progress
            if nargin ~= 0
                cost_new = calc_F(x_new);
                if cost_new <= cost_old 
                    stt = 'YES.';
                else 
                    stt = 'NO, check your code.';
                end
                fprintf('iter = %3d, cost = %d, cost decreases? %s\n', iter, cost_new, stt);
            else 
                if mod(iter, 5) == 0
                    fprintf('.');
                end
                if mod(iter, 10) == 0 
                   fprintf('%d', iter);
                end     
%                 if mod(iter, 100) == 0 
%                     fprintf('\n');
%                 end
            end        
        end 
        
    end
    X = x_new;
    %%
    if nargout == 2
        cost = calc_F(x_new);
    end
    %%
    if nargin == 0 
        mycost = calc_F(X)
        pause
    end
end 