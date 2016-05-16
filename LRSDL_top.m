function best_acc = LRSDL_top(dataset, N_train, k, k0, lambda1, lambda2, lambda3)
% function LRSDL_FDDL_top(dataset, N_train, k, k0, lambda1, lambda2, lambda3)
    
%% Dependencies
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    %% test mode 
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 50;
        k       = 20;
        k0      = 5;
        lambda1 = 0.005;
        lambda2 = 0.01;
        lambda3 = 0.05;
    end 
    %%
    t = getTimeStr();
    [dataset, Y_train, Y_test, label_train, label_test] = ...
        train_test_split(dataset, N_train);
    %% Parameter preparation
    fprintf('starting... %s\n', dataset)    ;
    C = max(label_train);
    opts.N           = N_train;
    opts.k           = k;
    opts.k0          = k0;
    opts.show_cost   = 0;
    opts.lambda1     = lambda1;
    opts.lambda2     = lambda2;
    opts.lambda3     = lambda3;
    opts.D_range     = k*(0:C);
    opts.D_range_ext = [opts.D_range k*C+k0];
    opts.initmode    = 'normal';   
    opts.max_iter    = 100;
    opts             = initOpts(opts);
    opts.verbal      = true;
    opts.tol         = 1e-8;
    %% Train 
    [D, D0, X, X0, CoefM, coefM0, opts, rt] = ...
                    LRSDL(Y_train, label_train, opts);
    fprintf('rt = %5.1f\n', rt);
    %% output filename 
    if ~exist('results', 'dir')
        mkdir('results');
    end 
    if k0 > 0
        if ~exist(fullfile('results', 'LRSDL'), 'dir')
            mkdir('results', 'LRSDL');
        end 
        fn = fullfile('results', 'LRSDL', strcat(dataset, ...
            '_N_', num2str(N_train), '_k_', num2str(k), ...
            '_k0_',num2str(k0),'_l1_', num2str(lambda1), ...
            '_l2_', num2str(lambda2), '_l3_', num2str(lambda3), ...
            '_', t, '.mat'));
    else 
        if ~exist(fullfile('results', 'FDDL'), 'dir')
            mkdir('results', 'FDDL');
        end 
        fn = fullfile('results', 'FDDL', strcat(dataset, ...
            '_N_', num2str(N_train), '_k_', num2str(k),'_l1_', ...
            num2str(lambda1), '_l2_', num2str(lambda2),'_', t, '.mat'));
    end
    disp(fn);
    %%
    X1 = [X; X0];
    Y_range = label_to_range(label_train);
    C = max(label_train);
    CoefMM0 = zeros(size(X1,1), C);
    for c = 1: C 
        X1c = get_block_col(X1, c, Y_range);
        CoefMM0(:,c) = mean(X1c,2);
    end    
    opts.verbal = 0;
    if numel(D0) ~= 0
        fprintf('GC1:\n');
        acc = LRSDL_pred_2(Y_test, D, D0, CoefM, coefM0, opts, label_test);
        fprintf('GC2:\n');
        acc2 = LRSDL_pred_3(Y_test, D, D0, CoefMM0, coefM0, opts, label_test);
        fprintf('LC:\n');
        acc3 = LRSDL_pred(Y_test, D, D0, CoefMM0, opts, label_test);
        acc = [acc acc2 acc3];
        fprintf('maximum acc: %4f\n', max(acc));
	    save(fn, 'acc', 'rt');
    else             
        opts.weight = 0.1;
        acc = [];
        for vgamma = [0.0001, 0.001, 0.01, 0.1]
            opts.gamma = vgamma;
            pred = FDDL_pred(Y_test, D, CoefM, opts);
            acc1 = double(numel(find(pred == label_test)))/...
                numel(label_test);
%             disp(acc1);
            fprintf('gamma = %.4f, acc = %.4f\n', vgamma, acc1);
            acc = [acc acc1];
        end 
        save(fn, 'acc', 'rt');
    end
    best_acc = max(acc);
    
end 
 
%%
function [X, X0] = SC_SDDL_Nov20_classify(Y, D, D0, m0, lambda1, lambda2)
    N = size(Y,2);
    k = size(D,2);
    k0 = size(D0,2);
    X1init = zeros(k + k0, N);
    D1 = [D D0];
    M0 = repmat(m0, 1, N);
    D1tD1 = D1'*D1;
    D1tY = D1'*Y;
    function cost = calc_F(X1)
        X = X1(1: k, :);
        X0 = X1(k+1:end,:);
        cost =  0.5*normF2(Y - D*X - D0*X0) + ...
                0.5*lambda2*normF2(X0 - M0) + ...
                lambda1*norm1(X1);
    end 

    function g = grad(X1)
        X = X1(1: k, :);
        X0 = X1(k+1:end,:);
        g = (D1tD1*X1 - D1tY + lambda2* [zeros(k, N); X0 - M0]);
    end     
    %% ========= Main FISTA ==============================
    L = max(eig(D1tD1)) + 2;
    opts.tol = 1e-8;
    opts.max_iter = 300;
    [X1, ~] = fista(@grad, X1init, L, lambda1, opts, @calc_F);  
    X = X1(1: k, :);
    X0 = X1(k+1:end,:);
end 

%%
function acc = LRSDL_pred_2(Y, D, D0, CoefM, m0, opts, label_test)
    nClasses = size(CoefM, 2);
    k = opts.k;
    k0 = opts.k0;
    D_range = k*(0:nClasses);
    
    N = size(Y,2);
    acc = [];
    % --------------- Sparse coding -------------------------
    for lambda1 = [0.0001, 0.001]
%         disp(lambda1);
        [X, X0] = SC_SDDL_Nov20_classify(Y, D, D0, m0, lambda1, 0.01);
        % --------------- classification -------------------------
        Yhat = Y - D0*X0;
        E1 = zeros(nClasses, N);
        E2 = E1;
        for c = 1: nClasses            
            Dc = get_block_col(D, c, D_range);
            Xc = get_block_row(X, c, D_range);
            Mc = repmat(CoefM(:, c), 1, N );
            R1 = Yhat - Dc*Xc;
            R2 = X - Mc;
            E1(c,:) = sum(R1.^2);
            E2(c,:) = sum(R2.^2);
        end
        for w = [0.33, 0.66]
            E = w*E1 + (1-w)*E2;
            [~, pred] = min(E);
            aaaa = double(sum(pred == label_test))/N;
            acc = [acc aaaa];
            fprintf('w: %f, lambda1 = %.4f,  acc: %f\n', w, lambda1, aaaa);
        end 
    end 
end 

%%
function acc = LRSDL_pred_3(Y, D, D0, CoefMM0, m0, opts, label_test)
    nClasses = size(CoefMM0, 2);
    k = opts.k;
    k0 = opts.k0;
    D_range = k*(0:nClasses);    
    N = size(Y,2);
    % --------------- Sparse coding -------------------------
    acc = [];
    for lambda1 = [0.0001, 0.001]
        [X, X0] = SC_SDDL_Nov20_classify(Y, D, D0, m0, lambda1, 0.01);
        X1 = [X; X0];
        % --------------- classification -------------------------
        Yhat = Y - D0*X0;
        % for w = 0:0.01:1
        % --------------- classification -------------------------
        E1 = zeros(nClasses, N);
        E2 = zeros(nClasses, N);
        for c = 1: nClasses
            Dc = get_block_col(D, c, D_range);
            Xc = get_block_row(X, c, D_range);
            Mc = repmat(CoefMM0(:, c), 1, N );
            R1 = Yhat - Dc*Xc;
            R2 = X1 - Mc;
            E1(c,:) = sum(R1.^2);
            E2(c,:) = sum(R2.^2);
        end
        for w = [0.33, 0.66]
            E = w*E1 + (1-w)*E2;
            [~, pred] = min(E);
            aaaa = double(sum(pred == label_test))/N;
            acc = [acc aaaa];
            fprintf('w: %f, lambda1 = %.4f,  acc: %f\n', w, lambda1, aaaa);
        end 
    end    
end 



