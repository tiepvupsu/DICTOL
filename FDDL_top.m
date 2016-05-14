function FDDL_top(dataset, N_train, k, lambda1, lambda2)
% function FDDL_top(dataset, N_train, k, lambda1, lambda2)
    %% Dependencies
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 50;
        k = 25;
        lambda1 = 0.001;
        lambda2 = 0.05;
    end 
    
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);
    %% Parameter preparation
    fprintf('starting... %s\n', dataset)    ;
    C                = max(label_train);
    k0               = 0;
    opts.N           = N_train;
    opts.k           = k;
    opts.k0          = 0;
    opts.show_cost   = 0;
    opts.lambda1     = lambda1;
    opts.lambda2     = lambda2;
    opts.lambda3     = 0;
    opts.D_range     = k*(0:C);
    opts.D_range_ext = [opts.D_range k*C+k0];
    opts.initmode    = 'normal';   
    opts.max_iter    = 100;
    opts             = initOpts(opts);
    opts.verbal      = true;
    opts.tol         = 1e-8;
    %% Train 
    [D, ~, X, X0, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
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
    %%
    X1 = [X; X0];
    Y_range = label_to_range(label_train);
    C = max(label_train);
    opts.verbal = 0;
    opts.weight = 0.1;
    acc = [];
    for vgamma = [0.0001, 0.001, 0.01, 0.1]
        opts.gamma = vgamma;
        pred = FDDL_pred(Y_test, D, CoefM, opts);
        acc1 = double(numel(find(pred == label_test)))/...
            numel(label_test);
        fprintf('gamma = %.4f, acc = %.4f\n', vgamma, acc1);
        acc = [acc acc1];
    end 
    save(fn, 'acc', 'rt');
end 
 