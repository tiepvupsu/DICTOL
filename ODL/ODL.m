function [D, X] = ODL(Y, k, lambda, opts, method)
% * Solving the following problem:
%  (D, X) = \arg\min_{D,X} 0.5||Y - DX||_F^2 + lambda||X||_1
% * Syntax: `[D, X] = ODL(Y, k, lambda, opts, sc_method)`
%   - INPUT: 
%     + `Y`: collection of samples.4/7/2016 7:35:39 PM
%     + `k`: number of atoms in the desired dictionary.
%     + `lambda`: norm 1 regularization parameter.
%     + `opts`: option.
%     + `sc_method`: sparse coding method used in the sparse coefficient update. Possible values:
%       * `'fista'`: using FISTA algorithm. See also [`fista`](#fista).
%       * `'spams'`: using SPAMS toolbox [[12]](#fn_spams). 
%   - OUTPUT:
%     + `D, X`: as in the problem.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/7/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0
		addpath(fullfile('..', 'utils'));
        addpath(fullfile('..', 'build_spams'));
        d      = 50; % data dimension
        N      = 100; % number of samples 
        k      = 50; % dictionary size 
        lambda = 0.1;
        Y      = normc(rand(d, N));        
        %
		opts.max_iter      = 500;
		opts.show_progress = 0;
		opts.check_grad    = false;  
		opts.tol           = 1e-8;  
		opts.verbose     = true;
		
	end 
	%%
	opts = initOpts(opts);
	%%
	%% ========= initial D ==============================
	D = PickDfromY(Y, [0, size(Y,2)], k);
    X = zeros(size(D,2), size(Y,2));
    if opts.verbose 
        fprintf('cost: %f', ODL_cost(Y, D, X, lambda));
    end 
    optsX = opts;
	optsX.max_iter = 200;
	optsX.tol      = 1e-8;
    optsD = opts;
	optsD.max_iter = 200;
	optsD.tol      = 1e-8;
	iter = 0;
	while iter < opts.max_iter
		iter = iter + 1;
		%% ========= sparse coding step ==============================
		X = lasso_fista(Y, D, X, lambda, optsX);
       	if opts.verbose 
			costX = ODL_cost(Y, D, X, lambda);
			fprintf('iter: %3d, costX = %5f\n', iter, costX)
		end 
		%% ========= dictionary update step ==============================
		F = X*X'; E = Y*X';
		D = ODL_updateD(D, E, F, optsD);
		if opts.verbose 
			costD = ODL_cost(Y, D, X, lambda);
			fprintf('iter: %3d, costD = %5f\n', iter, costD)
		end 
	end
	%%
	if nargin == 0
		pause;
	end
end
