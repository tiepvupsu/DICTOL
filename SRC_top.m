function best_acc = SRC_top(dataset, N_train, lambda)
% function SRC_top(dataset, N_trn, lambda)
% Description       : SRC 
%     INPUT: 
%       dataset: name of the dataset stored in 'data', excluding '.mat'
%       N_trn: number of training images per class 
%       lambda : regularization parameter lambda     
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    addpath('utils');
    addpath('SRC');
    addpath('build_spams');
    %% ========= Test mode ==============================
    if nargin == 0         
        dataset = 'myYaleB';
        N_train = 10;        
        lambda = 0.001;
    end 
    %%
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    %% output file 
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'SRC'), 'dir')
        mkdir('results', 'SRC');
    end 
    fn = fullfile('results', 'SRC', strcat(dataset, '_N_', num2str(N_train), ...
        '_l_', num2str(lambda), '_', t, '.mat'));
    %% main 
    range_train = label_to_range(label_train);
    range_test = label_to_range(label_test);
    acc = SRC_wrapper(Y_train, range_train, Y_test, range_test, lambda);   
    save(fn, 'acc');
    best_acc = acc;
end 
