function acc = LRSDL_pred_GC(Y, D, D0, CoefM, m0, opts, label_test)
    nClasses = size(CoefM, 2);
    k = opts.k;
    k0 = opts.k0;
    D_range = k*(0:nClasses);
    
    N = size(Y,2);
    acc = [];
    % --------------- Sparse coding -------------------------
    for lambda1 = [0.001]
        [X, X0] = local_sparse_coding(Y, D, D0, m0, lambda1, 0.01);
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
        for w = [ .5]
            E = w*E1 + (1-w)*E2;
            [~, pred] = min(E);
            aaaa = double(sum(pred == label_test))/N;
            acc = [acc aaaa];
            fprintf('w: %f, lambda1 = %.4f,  acc: %f\n', w, lambda1, aaaa);
        end 
    end 
end 
%%
function [X, X0] = local_sparse_coding(Y, D, D0, m0, lambda1, lambda2)
    N      = size(Y,2);
    k      = size(D,2);
    k0     = size(D0,2);
    X1init = zeros(k + k0, N);
    D1     = [D D0];
    M0     = repmat(m0, 1, N);
    D1tD1  = D1'*D1;
    D1tY   = D1'*Y;
    %% cost
    function cost = calc_F(X1)
        X = X1(1: k, :);
        X0 = X1(k+1:end,:);
        cost =  0.5*normF2(Y - D1*X1) + ...
                0.5*lambda2*normF2(X0 - M0) + ...
                lambda1*norm1(X1);
    end 
    %% grad
    function g = grad(X1)
        X  = X1(1: k, :);
        X0 = X1(k+1:end,:);
        g  = (D1tD1*X1 - D1tY + lambda2* [zeros(k, N); X0 - M0]);
    end     
    %% ========= Main FISTA ==============================
    L             = max(eig(D1tD1)) + 2;
    opts.tol      = 1e-8;
    opts.max_iter = 300;
    [X1, ~]       = fista(@grad, X1init, L, lambda1, opts, @calc_F);  
    X             = X1(1: k, :);
    X0            = X1(k+1:end,:);
end 

