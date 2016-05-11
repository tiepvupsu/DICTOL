function [D, X] = myDLCOPAR(Y, Y_range, opts)

 % cost(X,D) = normF2(Y - D*X) + sum_c((normF2(Yc - DcXcc - D0X0c) + sum_{i\neq c} normF2(Xic)) + lambda*norm1(X) + ...
 %    eta*sum_c sum_{i = c+1 -> C+1} normF2(Di'*Dc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ================== block: test module 12/19/2015 2:19:37 PM - Tiep Vu ==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%     function myDLCOPAR()
if nargin == 0
 profile off;
    profile on;
    C = 100;    N = 7;    d = 300;
    % C = 3;    N = 30;    d = 78;
    k = 7;
    k0 = 10;
    opts.k = k;
    opts.k0 = k0;
    opts.lambda = 0.0001;
    opts.eta = 0.01;
    opts.gamma = 0.01;
    
    opts.max_iter = 3;
    opts.show = false;
    Y = normc(rand(d, C*N));
    train_label = [];
    for c = 1: C
        train_label = [train_label c*ones(1, N)];
    end 
    D = normc(rand(d, C*k + k0));
    % D0 = normc(rand(d, k0));
    X = zeros(size(D,2), size(Y,2));
    % X0 = zeros(size(D0,2), size(Y,2));
    Y_range = N*(0:C);
    opts.D_range = [k*(0:C)];
    opts.max_iter = 5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ------------------end of block: test module 12/19/2015 2:19:37 PM - Tiep Vu ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    C = numel(Y_range) - 1;
    lambda = opts.lambda;
    eta = opts.eta;
    D_range = opts.D_range;
    % --------------- append D0 -------------------------
    D_range_ext = [D_range D_range(end)+opts.k0];
    k0 = opts.k0;
    opts.D_range_ext = D_range_ext;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: init ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin > 0
        d = size(Y,1);
        D = zeros(d, D_range_ext(end));
        X = zeros(D_range_ext(end), Y_range(end));
        parsinit.max_iter = 10;
        parsinit.show = false;
        for c = 1: C 
            Yc = get_block_col(Y, c, Y_range);
            [Dc, Xcc] = myODL(Yc, D_range(c+1) - D_range(c), lambda, parsinit);
            D(:, D_range(c)+1: D_range(c+1)) = Dc;
            X(D_range(c)+1: D_range(c+1), Y_range(c)+1: Y_range(c+1)) = Xcc;
        end 
        [DCp1, XCp1] = myODL(Y, k0, lambda, parsinit);
        D(:, D_range_ext(C+1)+1: end) = DCp1;
        X(D_range_ext(C+1)+1:end, :) = XCp1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: init ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    iter = 0;
    while iter < opts.max_iter 
        tic
        iter = iter + 1;

        %% ========= update X ==============================
        X = DLCOPAR_updateX(Y, Y_range, D, X, opts);
                
        fprintf('iter = %3d/%3d || costX = %5f\n', iter, opts.max_iter, DLCOPAR_cost(Y, Y_range, D, X, opts));        
        %% ========= update D ==============================
        D = DLCOPAR_updateD(Y, Y_range, D, X, opts); % and DCp1 
                
        fprintf('iter = %3d/%3d || costD = %5f\n', iter, opts.max_iter, DLCOPAR_cost(Y, Y_range, D, X, opts));        
        %% ========= update DCp1 ==============================
        
        % fprintf('iter = %3d || costDCp1 = %5f\n', iter, calccost(D,X));        
        %% ========= checkcost ==============================
        % cost = calccost(D, X);
        toc
    end 
end 




