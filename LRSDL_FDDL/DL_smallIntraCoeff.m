function [D, X] = DL_smallIntraCoeff(Y, k, lambda1, lambda2, opts)
% function [D, X] = DL_smallIntraCoeff(Y, k, lambda1, lambda2, opts)
% cost(D, X) = 0.5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildMean(X)) 
%               + lambda1*norm1(X);            
    if nargin == 0
    	d = 3000;
    	k = 10;
    	N = 10;
    	Y = normc(rand(d, N));
    	lambda1 = 0.001;
    	lambda2 = 0.5;
        opts.verbose = 0;
        opts.max_iter = 30;
    end 	
    opts = initOpts(opts);
    %%
	D = PickDfromY(Y, [0, size(Y,2)], k);
	X = zeros(size(D,2), size(Y,2));
    %%
	function cost = calc(D, X)
		cost = 0.5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildMean(X)) + lambda1*norm1(X);		
    end
    %%
	it = 0;
    cost_old = calc(D, X);    
    if opts.verbose
        fprintf('%f', cost_old);
        fprintf('.');
    end 
    %%
    optsX = opts;
    optsX.max_iter = 100;
    optsD = opts;
    optsD.max_iter = 30;
    while it < opts.max_iter
		it = it + 1;
		%% ========= update X ==============================
        % X = \arg\min_X 0.5*normF2(Y - D*X) + 
        %     .5*lambda2*normF2(X - buildMean(X)) + lambda1*norm1(X);
		X = myLassoWIntrasmall_fista(Y, D, lambda1, lambda2, X, opts);
		%% ========= update D ==============================
		E = Y*X';
		F = X*X';
		D = ODL_updateD(D, E, F, optsD);
		%% ========= For debugging ==============================
%         cost_new = calc(D, X);
        %if(abs(cost_new - cost_old) < 1e-6);
        %    break;
        % end
        % cost_old = cost_new;
    end 
    if opts.verbose
        fprintf('%f', cost_new);
    end 
    if nargin == 0
        X = [];
        D = [];
    end
end 

