function [D, X, rt] = DLCOPAR(Y, Y_range, opts)
    % syntax = [D, X] = DLCOPAR(Y, Y_range, opts)

 % cost(X,D) = normF2(Y - D*X) + sum_c((normF2(Yc - DcXcc - D0X0c) + sum_{i\neq c} normF2(Xic)) + lambda*norm1(X) + ...
     %    eta*sum_c sum_{i = c+1 -> C+1} normF2(Di'*Dc);
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
    
    C = numel(Y_range) - 1;
    lambda = opts.lambda;
    eta = opts.eta;
    D_range = opts.D_range;
    % --------------- append D0 -------------------------
    D_range_ext = [D_range D_range(end)+opts.k0];
    k0 = opts.k0;
    opts.D_range_ext = D_range_ext;
    
    %% ================== block: init ==========================    
    if nargin > 0
        d = size(Y,1);
        D = zeros(d, D_range_ext(end));
        X = zeros(D_range_ext(end), Y_range(end));
        parsinit.max_iter = 10;
        parsinit.show = false;
        for c = 1: C 
            fprintf(' class %3d: ', c);
            Yc = get_block_col(Y, c, Y_range);
            [Dc, Xcc] = myODL(Yc, D_range(c+1) - D_range(c), lambda, parsinit);
            D(:, D_range(c)+1: D_range(c+1)) = Dc;
            X(D_range(c)+1: D_range(c+1), Y_range(c)+1: Y_range(c+1)) = Xcc;
        end 
        [DCp1, XCp1] = myODL(Y, k0, lambda, parsinit);
        D(:, D_range_ext(C+1)+1: end) = DCp1;
        X(D_range_ext(C+1)+1:end, :) = XCp1;
    end
    %%
    iter = 0;
    optsX = opts;
    optsX.verbal = 0;
    optsD = opts;
    optsD.verbal = 0;
    tic
    while iter < opts.max_iter 
        iter = iter + 1;
        %% ========= update X ==============================
        if opts.verbal
            fprintf('Updating X...');
        end 
        X = DLCOPAR_updateX(Y, Y_range, D, X, optsX);
        t = toc;     
        if t > 20*3600
            break;
        end 
        if opts.verbal           
            costX = DLCOPAR_cost(Y, Y_range, D, X, opts);
            fprintf('iter = %3d || costX = %5f, running time: %f (s)\n', ...
                iter, costX, t);        
        end 
        %figure(1);
        %imagesc(X);
        %colormap jet
        %drawnow();
        %% ========= update D ==============================
        if opts.verbal
            fprintf('Updating D...');
        end 
        D = DLCOPAR_updateD(Y, Y_range, D, X, optsD); % and DCp1 
        if opts.verbal
            costD = DLCOPAR_cost(Y, Y_range, D, X, opts);        
            fprintf('iter = %3d || costD = %5f, running time: %f (s)\n', ...
                iter, costD, t);    
            if abs(costX - costD) < 1e-3
                break;
            end
        end 
        %figure(2);    
        %display_network(D);
        %drawnow();
        t = toc;       
        if t > 20*3600
            break;
        end
    end 
    rt = toc;
end 




