function DLCOPAR_top(dataset, N_train, k, k0, lambda, eta)
    myinit(); 
    if nargin == 0
        C = 5;
        d = 1000;
        N_train = 20;
        N_test = 10;
        k = 10;
        k0 = 15;
        D_range = k*(0:C);
        train_range = N_train*(0:C);
        test_range = N_test*(0:C);
        % D = normc(rand(d, D_range(end)));   
        Y_train = normc(rand(d, train_range(end)));   
        Y_test = normc(rand(d, test_range(end)));   

        opts.D_range = D_range;
        opts.lambda = 0.01;
        opts.eta = 0.05;
        label_test = RangeToLabel(test_range);  
    else 
        [Y_train, label_train, Y_test, label_test] = pickTrainTest_2(dataset, N_train);
    end

    opts.k = k;
    opts.k0 = k0;
    C = max(label_test);
    D_range = k*(0:C);
    opts.lambda = lambda;
    opts.eta = eta;
    % opts.gamma = 0.001; % for classification 
    opts.D_range = D_range;
    opts.D_range_ext = [D_range k*C+k0];

    train_range = label_to_range(label_train);
    opts.show = false;
    opts.max_iter = 100;        
    %% ========= Train ==============================
    [D, X] = myDLCOPAR(Y_train, train_range, opts);
    %% ========= test ==============================
    acc = [];
    fn = fullfile('results', 'DLCOPAR', strcat(dataset, '_N_', num2str(N_train), '_k_', num2str(k), '_k0_', num2str(k0), '_l_', num2str(lambda), '_e_', num2str(eta), '.mat'));
    for vgamma = [0.001, 0.005, 0.01, 0.1]
        opts.gamma = vgamma
        opts.classify_mode = 'LC';
        pred = DLCOPAR_classify(Y_test, D, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('LC mode, acc = %5f\n', acc(end));

        opts.classify_mode = 'GC';
        pred = DLCOPAR_classify(Y_test, D, opts);
        acc = [acc double(sum(pred == label_test))/numel(label_test)];
        fprintf('GC mode, acc = %5f\n', acc(end));

        save(fn, 'acc');

    end 
end 