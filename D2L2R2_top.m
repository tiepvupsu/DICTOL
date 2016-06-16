function best_acc = D2L2R2_top(dataset, N_train, k, lambda1, lambda2, alpha)
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
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;        
        lambda1 = 0.001;
        lambda2 = 0.01;
        alpha = 0.01;        
    end 
    %% Data preparation
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    %% main 
    [acc, rt] = D2L2R2_wrapper(Y_train, label_train, Y_test, label_test,...
                        k, lambda1, lambda2, alpha);
    %% save results 
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
    
    save(fn, 'acc', 'rt');
    best_acc = max(acc);
end 