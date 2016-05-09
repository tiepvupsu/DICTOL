function [X, iter] = FDLDL_updateX(Y, D, X, W, Y_range, opts)
    
    if nargin == 0 
%         profile off;
%         profile on;
        addpath('utils')
        C = 3; k = 30; d = 30; N = 30; opts.checkgrad = 0;
        C = 5; k = 5; d = 5; N = 5; opts.checkgrad = 1;
        opts.K = C*k;
        Y = normc(rand(d, N*C));
        D = normc(rand(d, k*C));
        X = zeros(size(D, 2), size(Y, 2));
        % W = zeros(k*C, C);
        
        W = rand(k*C, C);      


        Y_range = N*(0:C);
%         for c = 1: C 
%             W(k*(c-1) + 1: k*c,c) = ones(k, 1);
%         end 
        opts.lambda1 = 0.1;
        opts.lambda2 = 0.1;
        opts.lambda3 = 0.01;
        
        opts.showprogress = 1;
        opts.max_iter = 30;
        opts.tol = 1e-6;

    end 
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    X = X;
    K = opts.K;
    C = size(W, 2);
    % X = Xinit;
    %% TODO 
    function cost = calc_f(X) 
        cost = lambda2*normF2(X);
        % cost0 = 0;
        m = mean(X, 2);
        C = size(W, 2);
        for c = 1: C 
            Yc   = get_block_col(Y, c, Y_range);
            Xc   = get_block_col(X, c, Y_range);
            wc   = W(:, c);
            
            D_c  = myMult_A_diag_b(D, wc);
            cost = cost + .5*normF2(Yc - D_c*Xc);
            % cost = cost + beta2*wc'*(ones(numel(wc), 1) - wc);
            for i = 1: C 
                if i == c 
                    continue;
                else 
                    Xi = get_block_col(X, i, Y_range);
                    cost = cost + 0.5*normF2(D_c*Xi);
                end 
            end 
            mc   = mean(Xc, 2);
            Mc   = repmat(mc, 1, size(Yc, 2));
            M    = repmat(m, 1, size(Yc, 2));
            cost = cost + lambda2*(normF2(Xc - Mc) - normF2(Mc - M));
        end 
    end 
    %%
    function cost = calc_F(X)
        cost = calc_f(X) + lambda1*norm1(X);
    end 
    %% TODO 

    %% ========= Compute BigD ==============================    
    DtD = D'*D;

%     tmp = ones(K, K);
%     for c = 1: C 
%         wc = W(:,c);
%         tmp = tmp + wc*wc';
%     end 
    tmp = W*W';% + ones(K, K);
    BigD = DtD.*tmp + 4*lambda2*eye(K);
    % DtY = D'*Y;
    DitYi = zeros(size(X));

    
    for c = 1: C 
        range_c = Y_range(c)+1: Y_range(c+1);
        wc = W(:, c);
        DitYi(:, range_c) = myMult_A_diag_b(D, wc)'*Y(:, range_c);
    end 
    
    Dtmp = -DitYi;% - DtY;
    %%
    function gradX = grad(X)
        gradX = BigD*X + Dtmp +  buildM_2Mbar(X, Y_range, lambda2);
    end 
    %%
    function g = num_grad(X)
        ep = 1e-4;
        g = zeros(size(X));
        for i = 1: size(X, 1)
            for j = 1: size(X,2)
                Xp = X;
                Xm = X;
                Xp(i,j) = Xp(i,j) + ep; fp = calc_f(Xp);
                Xm(i,j) = Xm(i,j) - ep; fm = calc_f(Xm);
                g(i,j)  = (fp - fm) / (2*ep);
            end 
        end 
    end 
    %%
    % if isfield(opts, 'checkgrad') && opts.checkgrad %check grad 
    %     fprintf('checking gradient...\n');
    %     X = rand(size(X));
    %     diff_grad = (grad(X) - num_grad(X));        
    %     fprintf('...done, the difference should be small, diff = %5f\n', normF2(diff_grad)); 
    %     fprintf('Press any key to continue...\n');
    %     pause
    % end 
     %% ========= Main FISTA ==============================
    L = 10*max(eig(BigD)) + 10;
    optsX.max_iter = 300;
    optsX.tol = 1e-8;
    if nargin == 0
        cost_before = FDLDL_cost(Y, D, X, W, Y_range, opts);
        fprintf('Cost before : %5f\n', cost_before);
    end 
    [X, iter] = fista(@grad, X, L, opts.lambda1, optsX, @calc_F);   
    if nargin == 0
        cost_after = FDLDL_cost(Y, D, X, W, Y_range, opts);
        fprintf('Cost after  : %5f\n', cost_after);
    end 
    % Linv = 1/L;    
    % lambdaLiv = lambda1*Linv;
    % x_old = Xinit;
    % y_old = Xinit;
    % t_old = 1;
    % iter = 0;
    % if opts.showprogress
    %     fprintf('(max_iter = %d): 0.', opts.max_iter);
    % end 
    % if nargin == 0 
    %     cost_old =  calc_F(x_old);
    % end 
    % opts.max_iter = 300;
    % while  iter < opts.max_iter
    %     iter = iter + 1;
    %     x_new = shrinkage(y_old - Linv*grad(y_old), lambdaLiv);
    %     t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
    %     y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
    %     % check stop criteria
    %     e = norm1(x_new - x_old)/numel(x_new);
    %     if e < opts.tol
    %         break;
    %     end
    %     % update
    %     x_old = x_new;
    %     t_old = t_new;
    %     y_old = y_new;
    %     % show progress
    %     if opts.showprogress
    %         if nargin ~= 0
    %             cost_new = calc_F(x_new);
    %             if cost_new <= cost_old 
    %                 stt = 'YES.';
    %             else 
    %                 stt = 'NO, check your code.';
    %             end
    %             fprintf('iter = %3d, cost = %d, cost decreases? %s\n', iter, cost_new, stt);
    %         else 
    %             if mod(iter, 2) == 0
    %                 fprintf('.');
    %             end
    %             if mod(iter, 10) == 0 
    %                fprintf('%d', iter);
    %             end     
    %             cost_new = calc_F(x_new);
    %             if cost_new <= cost_old 
    %                 stt = 'YES.';
    %             else 
    %                 stt = 'NO, check your code.';
    %             end
    %             fprintf('iter = %3d, cost = %d, cost decreases? %s\n', iter, cost_new, stt);
    %             if mod(iter, 100) == 0 
    %                 fprintf('\n');
    %             end
    %         end        
    %     end 
        
    % end
    % X = x_new;
    
    if nargin == 0
%         profile viewer;        
        pause         
    end

end 

function g2 = buildM_2Mbar(X, Y_range, lambda2)
    g2 = zeros(size(X));
    C = numel(Y_range) - 1;
    m = mean(X, 2);
    lambda4 = 2*lambda2;
    for c = 1: C 
        Xc = get_block_col(X, c, Y_range);
        mc = mean(Xc, 2);
        g2(:, Y_range(c) + 1: Y_range(c+1)) = repmat(lambda4*(m - 2*mc), 1, Y_range(c+1) - Y_range(c));
    end 
end 