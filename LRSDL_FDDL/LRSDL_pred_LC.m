function acc = LRSDL_pred_LC(Y, D, D0, CoefMM0, opts, label_test)
    nClasses = size(CoefMM0, 2);
    N = size(Y,2);    
    acc = [];
    optsX.max_iter = 500;
    for lambda = [0.001, 0.005, 0.01, 0.05] 
        E = zeros(nClasses, size(Y,2));
        %% ========= LC ==============================
        for c = 1: nClasses
            Di = [get_block_col(D, c, opts.D_range) D0];
            [Xi, ~] = lasso_fista(Y, Di, [], lambda, optsX);
            R = Y - Di*Xi;
            E(c,:) = 0.5*sum(R.^2,1) + lambda*sum(abs(Xi),1);
        end 
        [~, pred] = min(E);
        aaaa = double(sum(pred == label_test))/N;
        acc = [acc aaaa];
        fprintf('lambda = %.4f, acc: %f\n', lambda, aaaa);        
    end 
end 