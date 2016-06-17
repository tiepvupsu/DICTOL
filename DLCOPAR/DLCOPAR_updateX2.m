function X = DLCOPAR_updateX(Y, Y_range, D, X, opts) 
% function X = DLCOPAR_updateX(Y, Y_range, D, X, opts)
% updating X in DLCOPAR. 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% ----------------------------------------------- 
    if nargin == 0
        addpath(fullfile('..', 'utils'));
        fprintf('Test most\n');
        d = 30;
        N = 30;
        k = 20;
        k0 = 0;
        C = 30 ;
        opts.k0 = k0;
        opts.lambda = 0.01;
        opts.eta = 0.1;
        Y = normc(rand(d, C*N));
        D = normc(rand(d, C*k + k0));
        X = zeros(size(D,2), size(Y,2));
        Y_range = N*(0:C);
        opts.D_range = k* (0:C);
        opts.D_range_ext = [opts.D_range opts.D_range(end)+k0];
        opts.verbal = true;
        opts.max_iter = 30;
        opts.check_grad = 0;
    end 
%     opts = initOpts(opts);
    %%
    C = numel(Y_range) - 1;
    DtD = D'*D;
    DtY = D'*Y;
    DCp1 = get_block_col(D, C+1, opts.D_range_ext);
    DCp1tDCp1 = DCp1'*DCp1;
    if opts.k0 > 0
        L = max(eig(DtD)) + max(eig(DCp1tDCp1));
    else 
        L = max(eig(DtD));
    end 
    optsX = opts;
    optsX.verbal = false;
    optsX.max_iter = 100;
    %%
    for c = 1: C
        Xc = get_block_col(X, c, Y_range);
        X(:, Y_range(c)+1: Y_range(c+1)) = ...
            DLCOPAR_updateXc2(DtD, DCp1tDCp1, DtY, Y_range, Xc, c, L, optsX);
        if opts.verbal
            costXc = DLCOPAR_cost(Y, Y_range, D, opts.D_range_ext, X, optsX);
            fprintf('class = %3d | costXc: %5f\n', c, costXc);
        end
    end
    %%
    if nargin == 0
        X = [];
    end
end 

