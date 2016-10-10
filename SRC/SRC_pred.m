function [pred, X] = SRC_pred(Y, D, D_range, opts)
% function [pred, X] = SRC_pred(Y, D, D_range, opts)
% * Classification based on SRC.
% * Syntax: `[pred, X] = SRC_pred(Y, D, D_range, opts)`
%   - INPUT:
%     + `Y`: test samples.
%     + `D`: the total dictionary. `D = [D_1, D_2, ..., D_C]` 
%		with `D_c` being the _c-th_ class-specific dictionary.
%     + `D_range`: range of class-specific dictionaries in `D`. 
%     + `opts`: options.
%       * `opts.lambda`: `lambda` for the Lasso problem.
%       * `opts.max_iter`: maximum iterations of fista algorithm. 
%       * others.
%   - OUTPUT:
%     + `pred`: predicted labels of test samples.
%     + `X`: solution of the lasso problem.
% Ref:
% 1. Wright, John, et al. "Robust face recognition via sparse representation."
%    Pattern Analysis and Machine Intelligence, IEEE Transactions on 31.2 (2009)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/6/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0
	    C = 3;    
	    N = 10;    
	    d = 100;
	    k = 7;
	    
	    Y = normc(rand(d, C*N));
	    D = normc(rand(d, C*k));
     	D_range = k* (0:C);
     	
	    opts.k = k;
	    opts.max_iter = 30;
	    opts.show_progress = false;
        opts.lambda = 0.01;
	end 
	if ~isfield(opts, 'lambda')
		opts.lambda = 0.01;
	end
	%%
    Y = normc(Y);
    D = normc(D);
	%% ========= Main code ==============================
	opts = initOpts(opts);
	%%
	C = numel(D_range) - 1;
	pred = zeros(1, size(Y,2));
	%%
	fprintf('sparse coding...');
	X = lasso_fista(Y, D, [], opts.lambda, opts);
	fprintf('done\n');
	%%
	E = zeros(C, size(Y,2));
	for i = 1:C
		Xi = get_block_row(X, i, D_range);
		Di = get_block_col(D, i, D_range);
		R = Y - Di*Xi;
		E(i,:) = sum(R.^2, 1);
	end
	[~, pred] = min(E);
	%%
	if nargin == 0
		return;
		% pause;
	end
end