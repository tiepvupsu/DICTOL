function DLSI_top(dataset, N_train, k, lambda, eta)


    myinit 
    [Y_train, label_train, Y_test, label_test] = pickTrainTest_2(dataset, N_train);

    C = max(label_train);
    D_range = k*(0:C);
    opts.lambda = lambda;
    opts.eta = eta;
    opts.D_range = D_range;

    train_range = LabelToRange(label_train);
    opts.show = 1;
    opts.max_iter = 50;        
    %% ========= Train ==============================
    [D, X] = DLSI(Y_train, train_range, opts);
    %% ========= test ==============================
    pred = DLSI_pred(Y_test, D, opts);
    acc = double(sum(pred == label_test))/numel(label_test)

    fn = fullfile('results', 'DLSI', strcat(dataset, '_N_', num2str(N_train), '_k_', ...
        num2str(k), '_l_', num2str(lambda), '_e_', num2str(eta), '.mat'));
    save(fn, 'acc');
end 
