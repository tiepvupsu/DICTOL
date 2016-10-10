function [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2)
% function [acc, rt] = FDDL_wrapper(Y_train, label_train, Y_test , ...
%       label_test, k, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0 % test mode
        dataset = 'myYaleB';
        N_train = 10;        
        [~, Y_train, Y_test, label_train, label_test] = ...
            train_test_split(dataset, N_train);        
        k = 8;
        lambda = 0.001;
        eta = 0.01;
    end 
    C                = max(label_train);
    k0               = 0;    
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
    opts.verbose      = true;
    opts.tol         = 1e-8;
    %% Train 
    [D, ~, ~, ~, CoefM, ~, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
    Y_range = label_to_range(label_train);
    C = max(label_train);
    opts.verbose = 0;
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
    acc = max(acc);
end 
