function [ pred, acc] = FDLDL_pred(Y, D, W, opts, train_label)
    lambda1 = opts.lambda1;

    C = opts.C;
    
    N = size(Y,2);
    optsX.max_iter = 300;
    % [X, ~] = myLasso_fista_2(Y, D, [], lambda1, optsX);
    cnt = 0;
    for lambda1 = [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.3]
        %% ========= global ==============================
        cnt = cnt + 1;
        mu = sum(W, 2);
        D1 = myMult_A_diag_b(D, mu);
        X = myLasso_spams(Y, D1, lambda1);
        % --------------- classification -------------------------
        E1 = zeros(C, N);
       
        for c = 1: C
            R = Y - myMult_A_diag_b(D, W(:, c))*X;
            E1(c, :) = sum(R.^2, 1);
        end
        [~, pred] = min(E1);
        if nargin == 5 % train_label included
            acc(cnt) = double(sum(pred == train_label))/numel(train_label);
            fprintf('lambda1 = %5f, acc GC = %5f\n', lambda1, acc(cnt));
        end

        %% ========= local ==============================
        cnt = cnt + 1;
        E2 = zeros(C, N);
        optsX.max_iter = 200;
        for c = 1: C
            Dc =  myMult_A_diag_b(D, W(:,c));
            X = myLasso_fista_2(Y, Dc, [], lambda1, optsX);
            R = Y - Dc*X;
            E2(c,:) = sum(R.^2, 1) + lambda1*sum(abs(X), 1);
        end 
        [~, pred] = min(E2);
        if nargin == 5 % train_label included
            acc(cnt) = double(sum(pred == train_label))/numel(train_label);
            fprintf('lambda1 = %5f, acc LC = %5f\n', lambda1, acc(cnt));
        end
    end
        

end 
