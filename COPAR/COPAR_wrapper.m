function [acc, rt] = COPAR_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda, eta)
% function [acc, rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda1, lambda2, lambda3)
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
        k0 = 5;
        lambda = 0.001;
        eta = 0.01;        
    end 
    opts.k = k;
    opts.k0 = k0;
    C                = max(label_test);
    D_range_ext = [k*(0:C), k*C + k0];
    opts.lambda      = lambda;
    opts.eta         = eta;    
    train_range      = label_to_range(label_train);
    opts.show        = false;
    opts.max_iter    = 10;        
    opts.verbose      = true;
    opts = initOpts(opts);
    %% ========= Train ==============================
    [D, X, rt] = COPAR(Y_train, train_range, opts);
    %% ========= test ==============================
    acc = [];
    for vgamma = [0.0001, 0.001, 0.005, 0.01]
        opts.gamma = vgamma;
        fprintf('gamma %5f\n', vgamma);
        opts.classify_mode = 'LC';
        pred = COPAR_pred(Y_test, D, D_range_ext, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('LC mode, acc = %5f\n', acc(end));

        opts.classify_mode = 'GC';
        pred = COPAR_pred(Y_test, D, D_range_ext, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('GC mode, acc = %5f\n', acc(end));        
    end
    acc = max(acc);
end 
