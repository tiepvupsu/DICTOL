function [D, X] = DLSI(Y, Y_range, opts)
% [D, X] = argmin_{D, X}(\sum 0.5*\|Y_i - D_i X_i\|_F^2) + \lambda*norm1(X) + 
%                   0.5*eta * \sum_{i \neq j} \|Di^T*Dj\|_F^2 
% Xi = arg\min 0.5*||Y_i - D_i X_i\|_F^2 + \lambda \|X_i\|
% Di = \arg\min \|Y_i - D_i X_i\|_F^2 + \eta \|D_{\i}^T D_i\|_F^2 
% 
% Ref: 
% Ramirez, Ignacio, Pablo Sprechmann, and Guillermo Sapiro. 
% "Classification and clustering via dictionary learning with structured 
% incoherence and shared features." _Computer Vision and Pattern Recognition 
% (CVPR), 2010 IEEE Conference on. IEEE_, 2010. 
%  (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5539964&tag=1)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/7/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
%         profile off;
%         profile on;
        addpath(fullfile('..','utils'));
        addpath(fullfile('..','ODL'));
        
        C = 100;    N = 7;    d = 3000;        
        C = 3;    N = 30;    d = 78;
        k = 7;
        opts.k0 = 10;
        opts.lambda = 0.01;
        opts.eta = 0.1;
        pars.max_iter = 3;
        pars.show = false;
        Y = normc(rand(d, C*N));
        train_label = [];
        for c = 1: C
            train_label = [train_label c*ones(1, N)];
        end         
        Y_range = label_to_range(train_label);
        opts.D_range = k* (0:C);
        opts.show = true;
        opts.max_iter = 10;
    end 
    %%
    D_range = opts.D_range;
    C = numel(D_range) - 1
    D = zeros(size(Y,1), D_range(end));
    % X = zeros(size(D,2), size(Y,2))

    %% ========= X should be stored in cell ==============================
    clear X;
    for i = 1: C 
        rows = D_range(i+1) - D_range(i);
        cols = Y_range(i+1) - Y_range(i);
        X{i} = zeros(rows, cols);
    end 
    lambda = opts.lambda;
    eta = opts.eta;
    %%
    %% ========= init ==============================
    optsinit = opts;
    optsinit.show = 0;
    fprintf('Cost = %5f\n', DLSI_cost(Y, Y_range, D, D_range, X, opts));
    fprintf('Initializing...');
    opts.max_iter = 50;
    for i = 1: C 
        if opts.show 
            fprintf('.');
        end
        fprintf('Class: %3d  ', i);
        Yi = get_block_col(Y,i,Y_range);
        [Di, X{i}] = ODL(Yi, D_range(i+1) - D_range(i), lambda, optsinit);
        D(:, D_range(i) + 1: D_range(i+1)) = Di;
    end 
    % if opts.show 
    fprintf('\n');
    fprintf('Cost = %5f\n', DLSI_cost(Y, Y_range, D, D_range, X, opts));
    %%
    iter = 0;
    optsX.max_iter = 300;
    optsX.show = false;
    optsD.show = false;
    optsD.max_iter = 200;
    %%
    while iter < opts.max_iter
        iter = iter + 1;
        %% ========= update X ==============================
        for i = 1: C 
            Yi = get_block_col(Y, i, Y_range);
            Di = get_block_col(D, i, D_range);
            X{i} = lasso_fista(Yi, Di, X{i}, lambda, optsX);
        end 
        if opts.show 
            fprintf('iter = %3d || costX = %5f\n', iter, DLSI_cost(Y, ...
                Y_range, D, D_range, X, opts));
        end 
        %% ========= update D ==============================
        for i = 1: C 
            D_comi = D;
            D_comi(D_range(i)+1: D_range(i+1)) = [];
            Di = D(:, D_range(i)+1: D_range(i+1));
            Yi = get_block_col(Y, i, Y_range);
            Di = DLSI_updateD(Yi, X{i}, Di, D_comi', eta, optsD);
            D(:, D_range(i)+1: D_range(i+1)) = Di;
        end 
        if opts.show 
            fprintf('iter = %3d || costD = %5f\n', iter, DLSI_cost(Y, ...
                Y_range, D, D_range, X, opts));
        end 
    end 
    %%
    if opts.show
        pause
    end 

end 