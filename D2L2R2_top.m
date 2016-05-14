function D2L2R2_top(dataset, N_train, k, lambda1, lambda2, alpha)
% function D2L2R2_top(dataset, N_train, k, lambda1, lambda2, alpha)
% * The top function of D2L2R2 
% * INPUT:
%   + `dataset`: name of the dataset stored in `.mat` file in `data` folder. 
%     Note that `dataset` is the file name of the `.mat`, excluding `.mat`.
%   + `N_train`: number of training samples in each class 
%   + `k`: number of atoms in EACH dictionary 
%   + `lambda1, lambda2, alpha`: regularization parameters.
% * To run an small example, type `D2L2r2-top` without input in MATLAB 
%           command window.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/13/2016 10:29:48 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath(genpath('utils'));  
    addpath('ODL');
    addpath('LRSDL_FDDL');
    addpath('D2L2R2');
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 40;
        k = 20;
        lambda1 = 0.001;
        lambda2 = 0.01;
        alpha = 0.01;
        % eta = 0.1;
    end 
    %% Data preparation
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'D2L2R2'), 'dir')
        mkdir('results', 'D2L2R2');
    end 
    fn = fullfile('results', 'D2L2R2', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_l1_', num2str(lambda1),...
        '_l2_', num2str(lambda2), '_a_', num2str(alpha), '_', t, '.mat'));
    disp(fn);
    %% parameter preparation 
    C             = max(label_train);    
    opts.k        = k;
    opts.lambda1  = lambda1;
    opts.lambda2  = lambda2;
    opts.alpha    = alpha;
    opts.max_iter = 100;        
    opts.show     = false;
    opts.showD    = false;
    opts.verbal   = true;
    opts          = initOpts(opts);
    %% ========= Train ==============================
    [D, D_range, X, CoefM, opts, rt] = D2L2R2(Y_train, label_train, opts);
    %% ========= test ==============================
    acc = [];
    opts.verbal = false;
    opts.max_iter = 300;
    for vgamma = [0.001, 0.005, 0.01, 0.1]
        opts.gamma = vgamma;
        opts.weight = 0.5;

        pred = D2L2R2_pred(Y_test, D, D_range, CoefM, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('gamma = %.4f, acc = %5f\n', vgamma, acc(end));
    end 
    save(fn, 'acc', 'rt');
end 