function X = myLassoWIntrasmall_fista(Y, D, lambda1, lambda2, Xinit, opts)
% function X = myLassoWIntrasmall_fista(Y, D, lambda1, lambda2, Xinit, pars)
% cost = .5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildM(X)) ...
%       + lambda1*norm1(X);
	if nargin == 0
        clc
        addpath(fullfile('..', 'utils'));
        d = 30;
        k = 7;
        N = 7;
        Y = normc(rand(d, N));
        D = normc(rand(d, k));
        
        lambda1 = 0.01;
        lambda2 = 1;

        Xinit = zeros(size(D,2), size(Y, 2));   
        opts.max_iter = 100;
        opts.verbose = true;
    end 
    %%
    opts = initOpts(opts);
	DtD = D'*D;
	DtY = D'*Y;
	I = eye(size(D,2));
	D1 = DtD + lambda2*I;
    %%
	function M = buildM(X)
		M = repmat(mean(X,2), 1, size(X,2));
    end
    %% cost   
    function cost = calc_f(X)
		cost = .5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildM(X));
    end 
    %%
	function cost = calc_F(X)
		cost = calc_f(X) +  lambda1*norm1(X);
    end 
    %% calculating gradient
	function grad = grad(X)
		grad = D1*X - DtY - lambda2*buildM(X);
    end 
    %% check grad
    if opts.check_grad
        check_grad(@calc_f, @grad, rand(size(X)));
    end     
    %% L 
    L = max(eig(D1)) + lambda2;
    X = fista(@grad, Xinit, L, lambda1, opts, @calc_F);
    %%
    if nargin == 0
        X = [];
    end 

end
