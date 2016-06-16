function best_acc = LCKSVD_top(dataset, N_train, k, ...
        sparsitythres, valpha, vbeta)
    % Syntax: LCKSVD_top(dataset, N_train, k, sparsitythres, valpha, vbeta)
    addpath(genpath('utils'));  % add K-SVD box
    addpath(genpath('LCKSVD'));  % add K-SVD box
    addpath('build_spams');
    addpath('utils');
    best_acc = zeros(1, 2);
    if nargin == 0 
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;        
        sparsitythres = 10;
        valpha = 0.01;
        vbeta = 0.01;
    end 
    %% get data 
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    %% main 
    [acc, rt] = LCKSVD_wrapper(Y_train, label_train, Y_test, label_test,...
                    k, sparsitythres, valpha, vbeta);
    %% save data
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if ~exist(fullfile('results', 'LCKSVD'), 'dir')
            mkdir('results', 'LCKSVD');
        end 
    fn1 = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), '_b_', ...
        num2str(vbeta), '_', getTimeStr(), '_1.mat'));
    disp(fn1);
    fn2 = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), ...
        '_b_', num2str(vbeta), '_', getTimeStr(), '_2.mat'));
    best_acc = acc;
    rrt = rt;
    acc = best_acc(1);
    rt = rrt(1);
    save(fn1, 'acc', 'rt');
    acc = best_acc(2);
    rt = rrt(2);
    save(fn2, 'acc', 'rt');
end 
