function [acc, rt] = D2L2R2_wrapper(Y_train, label_train, Y_test, label_test,...
                        k, lambda1, lambda2, alpha)
    if nargin == 0 %% test mode
        addpath(fullfile('..', 'utils'));  
        addpath(fullfile('..', 'ODL'));
        addpath(fullfile('..', 'LRSDL_FDDL'));
        addpath(fullfile('..', 'D2L2R2'));
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;        
        lambda1 = 0.001;
        lambda2 = 0.01;
        alpha = 0.01;
        % eta = 0.1;
        [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
            dataset, N_train);
    end 
    %% Data preparation
    C             = max(label_train);    
    opts.k        = k;
    opts.lambda1  = lambda1;
    opts.lambda2  = lambda2;
    opts.alpha    = alpha;
    opts.max_iter = 100;        
    opts.show     = false;
    opts.showD    = false;
    opts.verbose   = true;
    opts          = initOpts(opts);
    %% ========= Train ==============================
    [D, D_range, X, CoefM, opts, rt] = D2L2R2(Y_train, label_train, opts);
    %% ========= test ==============================
    acc = [];
    opts.verbose = false;
    opts.max_iter = 300;
    for vgamma = [0.001, 0.005, 0.01, 0.1]
        opts.gamma = vgamma;
        opts.weight = 0.5;

        pred = D2L2R2_pred(Y_test, D, D_range, CoefM, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('gamma = %.4f, acc = %5f\n', vgamma, acc(end));
    end 
    best_acc = max(acc);
end 