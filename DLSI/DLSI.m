function [D, X, rt] = DLSI(Y, Y_range, opts)
% * function `[D, X, rt] = DLSI(Y, Y_range, opts)`
% * The main DLSI algorithm 
% * INPUT: 
%   - `Y, Y_range`: training samples and their labels 
%   - `opts`: 
%     + `opts.lambda, opts.eta`: `lambda` and `eta` in the cost function 
%     + `opts.max_iter`: maximum iterations. 
% * OUTPUT:
%   - `rt`: total running time of the training process.   
% * `[D, X, rt] = argmin_{D, X}(\sum 0.5*\|Y_c - D_c X_c\|_F^2) + 
%           \lambda*norm1(X) + 0.5*eta * \sum_{i \neq c} \|D_i^T D_c\|_F^2`
% * updating `X`:
%   `Xi = arg\min 0.5*||Y_c - D_c X_c\|_F^2 + \lambda \|X_c\|`
% * updating `D`:
%   `Di = \arg\min \|Y_c - D_c X_c\|_F^2 + \eta \|D_{\c}^T D_c\|_F^2`
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/14/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0 % test mode
        clc
        addpath(fullfile('..', 'utils'));
        d = 300;
        C = 10;
        N = 10;
        k = 5;
        Y = normc(rand(d, N*C));
        Y_range = N*(0:C);
        D_range = k*(0:C);
        
        opts.D_range = D_range;
        opts.k = 5;
        opts.lambda = 0.001;
        opts.eta = 0.1;
        opts.verbose = false;
        opts = initOpts(opts);        
    end 
    %%
    D_range = opts.D_range;
    C = numel(Y_range) - 1;
    D = zeros(size(Y,1), D_range(end));
    %% ========= X should be stored in cell since different sizes =============
    clear X;
    for c = 1: C 
        n_rows = D_range(c+1) - D_range(c);
        n_cols = Y_range(c+1) - Y_range(c);
        X{c} = zeros(n_rows, n_cols);
    end 
    lambda = opts.lambda;
    eta = opts.eta;    
    %% ========= init ==============================    
    if opts.verbose
        fprintf('Cost = %.4f\n', DLSI_cost(Y, Y_range, D, D_range, X, opts));
        fprintf('Initializing...\n');
        fprintf('class:\n');
    end 
    optsinit = opts;
    optsinit.verbose = 0;
    optsinit.max_iter = 50;
    for c = 1: C        
        if opts.verbose
            fprintf('%3d  ', c);
            if mod(c, 10) == 0
                fprintf('\n');
            end 
        end 
        Yc = get_block_col(Y,c,Y_range);        
        [Dc, X{c}] = ODL(Yc, D_range(c+1) - D_range(c), lambda, optsinit, ...
            'fista');
        D(:, D_range(c) + 1: D_range(c+1)) = Dc;
    end 
    % if opts.show 
    if opts.verbose
        fprintf('\n');
        fprintf('Cost = %.4f\n', DLSI_cost(Y, Y_range, D, D_range, X, opts));
    end 
    %%
    iter = 0;
    optsX = opts;
    optsX.max_iter = 300;
    optsX.verbose = false;
    
    optsD = opts;
    optsD.verbose = false;
    optsD.max_iter = 200;
    %%
    costD = 0;
    costX = 0;
    tic;
    while iter < opts.max_iter
        iter = iter + 1;
        %% ========= update X ==============================
        for c = 1: C 
            Yc = get_block_col(Y, c, Y_range);
            Dc = get_block_col(D, c, D_range);
            X{c} = lasso_fista(Yc, Dc, X{c}, lambda, optsX);
        end 
        if opts.verbose
            costX = DLSI_cost(Y, Y_range, D, D_range, X, opts);
            fprintf('iter = %3d/%3d | costX = %.4f\n', iter, ...
                opts.max_iter, costX);
        end 
        %% ========= update D ==============================
        for c = 1: C 
            D_comi = D;
            D_comi(:, D_range(c)+1: D_range(c+1)) = [];
            % Dc = D(:, D_range(c)+1: D_range(c+1));
            Dc = get_block_col(D, c, D_range);
            Yc = get_block_col(Y, c, Y_range);
%             Di = DLSI_updateD(Yi, X{i}, Di, D_comi', eta, optsD);
            E = Yc*X{c}';
            F = X{c}*X{c}';
            A = D_comi';
            Dc = DLSI_updateD(Dc, E, F, A, eta, optsD);
            D(:, D_range(c)+1: D_range(c+1)) = Dc;
        end 
        %%
        t0 = toc;
        if opts.verbose 
            costD = DLSI_cost(Y, Y_range, D, D_range, X, opts);
            fprintf('                 costD = %.4f', costD);
            t = t0*(opts.max_iter - iter)/iter;
            time_estimate(t);
        end 
        if t0 > 20*3600 % > 20h 
            break;
        end 
    end 
    rt = toc;
    if nargin == 0
        D = [];
        X = [];
    end 
end 
