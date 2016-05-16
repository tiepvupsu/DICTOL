function best_acc = DLSI_top(dataset, N_train, k, lambda, eta)
% * function DLSI_top(dataset, N_train, k, lambda, eta)
% * The top function of DLSI 
% * INPUT:
%   + `dataset`: name of the dataset stored in `.mat` file in `data` folder. 
%     Note that `dataset` is the file name of the `.mat`, excluding `.mat`.
%   + `N_train`: number of training samples in each class 
%   + `k`: number of atoms in EACH dictionary 
%   + `lambda, eta`: regularization parameters.
% * To run an small example, type `DLSI_top` without input in MATLAB command window. 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath(genpath('utils'));  
    addpath(genpath('DLSI'));      
    addpath('ODL')
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 50;
        k = 25;
        lambda = 0.001;
        eta = 0.1;
    end 
    %%
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'DLSI'), 'dir')
        mkdir('results', 'DLSI');
    end 
    fn = fullfile('results', 'DLSI', strcat(dataset, '_N_', num2str(N_train), ...
        '_k_', num2str(k), '_l_', num2str(lambda), '_e_', num2str(eta), '_', ...
        t, '.mat'));
    disp(fn);
    
    C              = max(label_train);
    D_range        = k*(0:C);
    opts.lambda    = lambda;
    opts.eta       = eta;
    opts.D_range   = D_range;
    opts.show_cost = 0;
    train_range    = label_to_range(label_train);
    opts.show      = 0;
    opts.verbal    = true;
    opts.max_iter  = 100;        
    %% ========= Train ==============================
    [D, X, rt]         = DLSI(Y_train, train_range, opts);
    %% ========= test ==============================
    opts.verbal    = false;
    pred           = DLSI_pred(Y_test, D, opts);
    acc            = double(sum(pred == label_test))/numel(label_test);
    disp(['acc = ', num2str(acc)]);
    disp(fn);    
    if strcmp(dataset, 'mySynthetic')
        save(fn, 'D', 'X', 'k', 'acc', 'rt');
    else 
        save(fn, 'acc', 'rt');
    end 
    best_acc = acc;
end 

