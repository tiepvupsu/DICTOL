function [acc, rt] = DLSI_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, lambda, eta)
% function acc = SRC_wrapper(Y_train, range_train, Y_test , range_test, lambda)
% Description       : SRC 
%     INPUT: 
%       dataset: name of the dataset stored in 'data', excluding '.mat'
%       N_trn: number of training images per class 
%       lambda : regularization parameter lambda     
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
    C              = max(label_train);
    D_range        = k*(0:C);
    opts.lambda    = lambda;
    opts.eta       = eta;
    opts.D_range   = D_range;
    opts.show_cost = 0;
    train_range    = label_to_range(label_train);
    opts.show      = 0;
    opts.verbose    = false;
    opts.max_iter  = 100;        
    %% ========= Train ==============================
    [D, ~, rt]         = DLSI(Y_train, train_range, opts);
    %% ========= test ==============================
    opts.verbose    = false;
    pred           = DLSI_pred(Y_test, D, opts);
    acc            = double(sum(pred == label_test))/numel(label_test);
end 
