function best_acc = COPAR_top(dataset, N_train, k, k0, lambda, eta)
% * function `COPAR_top(dataset, N_train, k, k0, lambda, eta)`
% * The top function of COPAR
% * INPUT:
%   + `dataset`: name of the dataset stored in `.mat` file in `data` folder. 
%     Note that `dataset` is the file name of the `.mat`, excluding `.mat`.
%   + `N_train`: number of training samples in each class 
%   + `k`: number of bases in EACH PARTICULAR dictionary 
%   + `k0`: number of bases in the COMMON dictionary
%   + `lambda, eta`: regularization parameters.
% * To run an small example, type `COPAR_top` without input in 
%     MATLAB command window.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath(genpath('utils'));  
    addpath(genpath('DLSI'));      
    addpath(genpath('COPAR'));      
    addpath('ODL')
    %% test mode 
    if nargin == 0 
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;
        k0 = 5;
        lambda = 0.001;
        eta = 0.01;
    end 
    %% get data 
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = ...
        train_test_split(dataset, N_train);
    %% main 
    [acc, rt] = COPAR_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda, eta);
    disp(rt);
    %% output filename 
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'COPAR'), 'dir')
        mkdir('results', 'COPAR');
    end 
    fn = fullfile('results', 'COPAR', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_k0_', num2str(k0), ...
        '_l_', num2str(lambda), '_e_', num2str(eta), '_', t, '.mat'));
    disp(fn);
    save(fn, 'acc', 'rt'); 
    best_acc = max(acc);
end 
