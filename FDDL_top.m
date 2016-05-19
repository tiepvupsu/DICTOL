function best_acc = FDDL_top(dataset, N_train, k, lambda1, lambda2)
% function FDDL_top(dataset, N_train, k, lambda1, lambda2)
    %% Dependencies
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 350;
        k = 25;
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;
        k0 = 5;
        lambda1 = 0.001;
        lambda2 = 0.05;
    end 
    %% get data 
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    %% main 
    [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2)
    fprintf('rt = %5.1f\n', rt);
    %% output filename 
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'FDDL'), 'dir')
        mkdir('results', 'FDDL');
    end 
    fn = fullfile('results', 'FDDL', strcat(dataset, ...
        '_N_', num2str(N_train), '_k_', num2str(k),'_l1_', ...
        num2str(lambda1), '_l2_', num2str(lambda2),'_', t, '.mat'));
    disp(fn);
    save(fn, 'acc', 'rt');
    best_acc = max(acc);
end 
 