function acc = SRC_wrapper(Y_train, range_train, Y_test , range_test, lambda)
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
        lambda = 0.001;
        [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
            dataset, N_train);
        range_train = label_to_range(label_train);
        range_test = label_to_range(label_test);
    end 
    %%
    opts.lambda   = lambda;
    opts.verbose   = false;
    opts.max_iter = 300;
    %%
    pred        = SRC_pred(Y_test, Y_train, range_train, opts);
    label_test = range_to_label(range_test);
    acc = calc_acc(pred, label_test);    
end 
