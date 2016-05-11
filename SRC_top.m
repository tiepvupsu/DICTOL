function SRC_top(dataset, N_trn, lambda)
    %% ================== File info ==========================
    % Author            : Tiep Vu (http://www.personal.psu.edu/thv102/)
    % Time created      : Thu Mar 10 10:07:57 2016
    % Last modified     : 
    % Description       : SRC 
    %     INPUT: 
    %       dataset: name of the dataset stored in 'data', excluding '.mat'
    %       N_trn: number of training images per class 
    %       lambda : regularization parameter lambda     
    %     OUTPUT: 
    %
    %% ================== end File info ==========================
    addpath('utils');
    addpath('basic_funcs');
    addpath('SRC');
    addpath('build_spams');

    %% ========= Test mode ==============================
    if nargin == 0
        dataset = 'myYaleB';
        N_trn = 15;
        lambda = 0.001;
    end 

    [Y_trn, label_trn, Y_tst, label_tst] = picktrntst_wrapper(dataset, N_trn);

    opts.lambda  = lambda ;
    opts.show    = 1;
    opts.max_iter= 300;
    
    trn_range = label_to_range(label_trn);
    pred  = SRC_pred(Y_tst, Y_trn, trn_range, opts);
    acc = double(sum(pred == label_tst))/numel(pred)

%     filename = fullfile('results','SRC', strcat('Normc_', dataset, '_N_', num2str(N_train), '_lambda_', num2str(lambda), '.mat'))
%     save(filename, 'acc');
end 
