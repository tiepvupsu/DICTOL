function best_acc = LRSDL_top(dataset, N_train, k, k0, lambda1, lambda2, lambda3)
% function LRSDL_FDDL_top(dataset, N_train, k, k0, lambda1, lambda2, lambda3)
    
%% Dependencies
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    %% test mode 
    if nargin == 0 
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;
        k0 = 5;
        lambda1 = 0.001;
        lambda2 = 0.01;
        lambda3 = 0.02;
    end 
    if ~exist('lambda3', 'var')
        lambda3 = 0;
    end 
    %%
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = ...
        train_test_split(dataset, N_train);
    %% Parameter preparation
    fprintf('starting... %s\n', dataset)    ;
    [acc, rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda1, lambda2, lambda3);
    fprintf('rt = %5.1f\n', rt);
    %% output filename 
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if k0 > 0
        if ~exist(fullfile('results', 'LRSDL'), 'dir')
            mkdir('results', 'LRSDL');
        end 
        fn = fullfile('results', 'LRSDL', strcat(dataset, ...
            '_N_', num2str(N_train), '_k_', num2str(k), ...
            '_k0_',num2str(k0),'_l1_', num2str(lambda1), ...
            '_l2_', num2str(lambda2), '_l3_', num2str(lambda3), ...
            '_', t, '.mat'));
    else 
        if ~exist(fullfile('results', 'FDDL'), 'dir')
            mkdir('results', 'FDDL');
        end 
        fn = fullfile('results', 'FDDL', strcat(dataset, ...
            '_N_', num2str(N_train), '_k_', num2str(k),'_l1_', ...
            num2str(lambda1), '_l2_', num2str(lambda2),'_', t, '.mat'));
    end
    disp(fn);
    save(fn, 'acc', 'rt');
    best_acc = max(acc);    
end 
 
  




