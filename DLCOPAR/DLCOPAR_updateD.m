function D = DLCOPAR_updateD(Y, Y_range, D, X, opts); % and DCp1 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: test ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0 
        C = 10;    N = 7;    d = 300;      
        k = 7;
        opts.k0 = 10;
        opts.lambda = 0.01;
        opts.eta = 0.1;
        pars.max_iter = 3;
        pars.show = false;
        Y = normc(rand(d, C*N));        
        Y_range = N * (0:C);
        opts.D_range = k* (0:C);
        opts.D_range_ext = [opts.D_range opts.D_range(end)+opts.k0];
        D = normc(rand(d, opts.D_range_ext(end)));
        opts.show = true;
        opts.max_iter = 10;
        X = 0.01*rand(size(D,2), size(Y,2));
    end         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: test ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    C = numel(Y_range) - 1;
    D_range_ext = opts.D_range_ext;
    for ii = 1: 1
    DCp1 = get_block_col(D, C+1, D_range_ext);
    optsD.max_iter = 100;
    optsD.show = false;    
    optsD.DLSI_YXt = true;
    optsD.lambda = 2*opts.eta;
%     fprintf('Updating D, cost = %5f\n', DLCOPAR_cost(Y, Y_range, D, X, opts));     
    Yhat = zeros(size(Y));
    for c = 1: C 
        %Dc = arg\min_Dc \|Ychat - Dc*Xcc\|_F^2 + \|Ycbar - Dc*Xcc\| + 2*eta\|A*Dc\|_F^2
        % = \arg\min_Dc \| [Ychat Ycbar] - Dc*[Xcc Xcc]\|_F^2 + 2*eta\|A*Dc\|_F^2
        % and solved using DLSI_updateD

        Dc_range = D_range_ext(c)+1: D_range_ext(c+1);
        Yc_range = Y_range(c)+1: Y_range(c+1);
        Yc = Y(:, Yc_range);
        Dc = D(:, Dc_range);
        Xc = X(:, Yc_range);
        Xcc = get_block_row(Xc, c, D_range_ext);
        XCp1c = get_block_row(Xc, C+1, D_range_ext);
%         Ytmp = Yc - Dc*Xcc;
        Ychat = Yc - D*Xc + Dc*Xcc;
        Ycbar = Yc - DCp1*XCp1c;
        YXt = (Ychat + Ycbar)*Xcc';
        XXt = 2*Xcc*Xcc';
        A = D;
        A(:,Dc_range) = []; 
        D(:, Dc_range) = DLSI_updateD(YXt, XXt, Dc, A', optsD);
%         fprintf('class = %3d| cost = %5f\n', c, DLCOPAR_cost(Y, Y_range, D, X, opts));
        
        
        Yhat(:, Yc_range) = Yc - D(:, Dc_range)*Xcc;    
    end 
    %%
    
%     for c = 1: C
%        Yc = 
%        Yhat(:, Yc_range) = Yc - Dc*Xcc;    
%     end

    %% ========= DCp1 ==============================
    XCp1 = X(D_range_ext(C+1) + 1 : D_range_ext(C+2), :);
    Ybar = Y - D(:, 1: D_range_ext(end-1))*X(1: D_range_ext(end-1), :);
    YXt = (Ybar + Yhat)*XCp1';
    XXt = 2*XCp1*XCp1';
    A = D(:, 1: D_range_ext(C+1));
    DCp1_range = D_range_ext(C+1) + 1: D_range_ext(C+2);
    D(:, DCp1_range) = DLSI_updateD(YXt, XXt, D(:, DCp1_range), A', optsD);         
%     fprintf('class = C+1| cost = %5f\n', DLCOPAR_cost(Y, Y_range, D, X, opts));
    end
end 