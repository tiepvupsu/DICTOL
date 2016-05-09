function pred = DLSI_classify(Y, D, opts)
    % j = \arg\min_j R(y, Dj) with R(y,D) = 0.5*\|y - Dx\|_2^2 + lambda*\|x\|_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ================== block: 12/18/2015 9:23:18 PM test module ==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DLSI_classify()
%     d = 50;
%     N = 200;
%     k = 100;
%     C = 10;
%     opts.D_range =  k*(0:C);
%     D = normc(rand(d, opts.D_range(end)));
%     opts.lambda = 0.01;
%     Y = normc(rand(d, N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ------------------end of block: 12/18/2015 9:23:18 PM test module ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    D_range = opts.D_range;
    C = numel(opts.D_range)-1;
    E = zeros(C, size(Y,2));
    pars.show = false;
    pars.max_iter = 300;
    for c = 1: C 
        Dc = get_block_col(D, c, D_range);
        Xc = myLasso_fista_2(Y, Dc, zeros(size(Dc,2), size(Y,2)), opts.lambda, pars);
        R1 = Y - Dc*Xc;
        E(c,:) = 0.5* sum(R1.^2, 1) + opts.lambda*sum(abs(Xc),1);
    end 
    [~, pred] = min(E);
end 
