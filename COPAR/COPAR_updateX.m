function X = COPAR_updateX(Y, Y_range, D, X, opts) 
% function X = COPAR_updateX(Y, Y_range, D, X, opts)
% updating X in COPAR. 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% ----------------------------------------------- 
    if nargin == 0
        addpath(fullfile('..', 'utils'));
        fprintf('Test most\n');
        d = 3000;
        N = 30;
        k = 20;
        k0 = 50;
        C = 100 ;
        opts.k0 = k0;
        opts.lambda = 0.01;
        opts.eta = 0.1;
        Y = normc(rand(d, C*N));
        D = normc(rand(d, C*k + k0));
        X = zeros(size(D,2), size(Y,2));
        Y_range = N*(0:C);
        opts.D_range = k* (0:C);
        opts.D_range_ext = [opts.D_range opts.D_range(end)+k0];
        opts.verbose = false;
        opts.max_iter = 30;
    end 
    %%
    C = numel(Y_range) - 1;
    DtD = D'*D;
    DtY = D'*Y;
    L = 2*max(eig(DtD)) + 10;    
    optsX = opts;
    optsX.max_iter = 100;
    %%
    for c = 1: C
        Xc = get_block_col(X, c, Y_range);
        X(:, Y_range(c)+1: Y_range(c+1)) = ...
            COPAR_updateXc(DtD, DtY, Y_range, Xc, c, L, optsX);
        if opts.verbose
            costXc = COPAR_cost(Y, Y_range, D, X, opts);
            fprintf('class = %3d | costXc: %5f\n', c, costXc);
        end
    end
    %%
    if nargin == 0
        X = [];
    end
end 

