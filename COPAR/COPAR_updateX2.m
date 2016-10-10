function X = COPAR_updateX2(Y, Y_range, D, X, opts) 
% function X = COPAR_updateX(Y, Y_range, D, X, opts)
% updating X in COPAR. 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% ----------------------------------------------- 
    if nargin == 0
        addpath(fullfile('..', 'utils'));
        fprintf('Test most\n');
        d = 30;
        N = 7;
        k = 7;
        k0 = 10;
        C = 20 ;
   
%         Y = normc(rand(d, C*N));
%         D = normc(rand(d, C*k + k0));
%         X = rand(size(D,2), size(Y,2));
%         Y_range = N*(0:C);
%         
%         save(fullfile('..', 'dictol_python', 'data', 'tmp5.mat'), ...
%             'Y', 'D', 'X', 'Y_range');
%         
        load(fullfile('..', 'dictol_python', 'data', 'tmp5.mat'), ...
            'Y', 'D', 'X', 'Y_range');
        
        opts.k0 = k0;
        opts.lambda = 0.01;
        opts.eta = 0.1;
        opts.D_range = k* (0:C);
        opts.D_range_ext = [opts.D_range opts.D_range(end)+k0];
        opts.verbose = true;
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
    optsX.verbose = false;
    optsX.max_iter = 100;
    %%
    for c = 1: C
        Xc = get_block_col(X, c, Y_range);
        X(:, Y_range(c)+1: Y_range(c+1)) = ...
            COPAR_updateXc2(DtD, DCp1tDCp1, DtY, Y_range, Xc, c, L, optsX);
        if opts.verbose
            costXc = COPAR_cost(Y, Y_range, D, opts.D_range_ext, X, optsX);
            fprintf('class = %3d | costXc: %5f\n', c, costXc);
        end
    end
    %%
    if nargin == 0
        X = [];
    end
end 

