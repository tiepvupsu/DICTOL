function pred = SRC_pred_shared(Y, D, D_range, opts)
% function SRC_pred_shared(Y, D, D_range, lambda, opts)
% SRC with presence of a shared class (ground, for example)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/19/2016 3:40:14 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 % test mode 
		addpath(fullfile('..', 'utils'));
		C = 3;
		N = 10;
		d = 100;
		k = 7;
		k0 = 5;

		Y = normc(rand(d, C*N));
	    D = normc(rand(d, C*k + k0));
     	D_range = [k*(0:C), k*C + k0];
     	
	    opts.k = k;
	    opts.max_iter = 30;
	    opts.show_progress = false;
        opts.lambda = 0.01;
    end 
    %%
    Y = normc(Y);
    D = normc(D);
    opts = initOpts(opts);
    if ~isfield(opts, 'norm')
        opts.norm = 1;
    end 
	%%
	C = numel(D_range) - 2;
	pred = zeros(1, size(Y,2));
	%%
	% fprintf('sparse coding...');
    if opts.norm == 1 
	   X = lasso_fista(Y, D, [], opts.lambda, opts);
    else % norm 0 
        X = omp(D'*Y, D'*D, opts.L);
    end
	% fprintf('done\n');
	%%
	E = zeros(C, size(Y,2));
	Ybar = Y - get_block_col(D, C+1, D_range)*get_block_row(X, C+1, D_range);
	for i = 1:C
		Xi = get_block_row(X, i, D_range);
		Di = get_block_col(D, i, D_range);
		R = Ybar - Di*Xi;
		E(i,:) = sum(R.^2, 1);
	end
	[~, pred] = min(E);
	%%
	if nargin == 0
		return;
		% pause;
	end
end
