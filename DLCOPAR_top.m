function DLCOPAR_top(dataset, N_train, k, k0, lambda, eta)
% * function `DLCOPAR_top(dataset, N_train, k, k0, lambda, eta)`
% * The top function of DLCOPAR
% * INPUT:
%   + `dataset`: name of the dataset stored in `.mat` file in `data` folder. 
%     Note that `dataset` is the file name of the `.mat`, excluding `.mat`.
%   + `N_train`: number of training samples in each class 
%   + `k`: number of bases in EACH PARTICULAR dictionary 
%   + `k0`: number of bases in the COMMON dictionary
%   + `lambda, eta`: regularization parameters.
% * To run an small example, type `DLCOPAR_top` without input in 
%     MATLAB command window.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath(genpath('utils'));  
    addpath(genpath('DLSI'));      
    addpath(genpath('DLCOPAR'));      
    addpath('ODL')
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 50;
        k = 25;
        k0 = 10;
        lambda = 0.001;
        eta = 0.01;
    end 
    %%
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = ...
        train_test_split(dataset, N_train);
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'DLCOPAR'), 'dir')
        mkdir('results', 'DLCOPAR');
    end 
    fn = fullfile('results', 'DLCOPAR', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_k0_', num2str(k0), ...
        '_l_', num2str(lambda), '_e_', num2str(eta), '_', t, '.mat'));
    %% options for DLCOPAR
    opts.k = k;
    opts.k0 = k0;
    C                = max(label_test);
    D_range_ext = [k*(0:C), k*C + k0];
    opts.lambda      = lambda;
    opts.eta         = eta;    
    train_range      = label_to_range(label_train);
    opts.show        = false;
    opts.max_iter    = 100;        
    opts.verbal      = true;
    opts = initOpts(opts);
    %% ========= Train ==============================
    [D, X, rt] = DLCOPAR(Y_train, train_range, opts);
    %% ========= test ==============================
    acc = [];
    for vgamma = [0.0001, 0.001, 0.005, 0.01]
        opts.gamma = vgamma;
        fprintf('gamma %5f\n', vgamma);
        opts.classify_mode = 'LC';
        pred = DLCOPAR_pred(Y_test, D, D_range_ext, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('LC mode, acc = %5f\n', acc(end));

        opts.classify_mode = 'GC';
        pred = DLCOPAR_pred(Y_test, D, D_range_ext, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('GC mode, acc = %5f\n', acc(end));
        if strcmp(dataset, 'mySynthetic')
            save(fn, 'D', 'X',  'opts', 'acc', 'rt');
        else 
            save(fn, 'acc', 'rt');
        end
    end
    disp(fn);
    save(fn, 'acc', 'rt'); 
end 
