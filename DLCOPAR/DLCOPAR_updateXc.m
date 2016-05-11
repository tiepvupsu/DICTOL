function Xc = DLCOPAR_updateXc(DtD, DtY,  Y_range, Xc, c, L, opts) %page 189-190 DLCOPAR
    % see DLCOPAR paper: http://www.cs.zju.edu.cn/people/wangdh/papers/draft_ECCV12_particularity.pdf
    % cost = normF2(Yc - D*Xc) + normF2(Yc - DcXcc - DCp1*XCp1c) + sum_{i \neq c, 1 \leq i \leq C} normF2(Xic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ================== block: test module 12/19/2015 4:03:53 PM =========================%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DLCOPAR_updateX()
% nargin
if nargin == 0
    fprintf('Test most\n');
    % profile off;
    % profile on;
    C = 100;    N = 10;    d = 30;
    
    % C = 3;    N = 30;    d = 78;
    k = 10;
    k0 = 10;    
    opts.k0 = k0;
    opts.lambda = 0.01;
    opts.eta = 0.1;
    opts.max_iter = 300;
    opts.show = true;
    Y = normc(rand(d, C*N));
    D = normc(rand(d, C*k + k0));
%     train_label = [];
%     for c = 1: C
%         train_label = [train_label c*ones(1, N)];
%     end         
    Y_range = N*(0:C);
    opts.D_range = k* (0:C);
    opts.D_range_ext = [opts.D_range opts.D_range(end)+k0];
    DtD = D'*D;
    DtY = D'*Y;
    c = 2;
    Yc = get_block_col(Y, c, Y_range);
    Xc = zeros(size(D,2), size(Yc,2));
    L = 2*max(eig(DtD))+10;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ------------------end of block: test module 12/19/2015 4:03:53 PM -------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% tic 
    C = numel(opts.D_range_ext) - 2;

    % d = size(Y,1);
    D_range_ext = opts.D_range_ext;  
    
    lambda = opts.lambda/2;
    function cost = calc_F(Xc)
        cost = normF2(Yc - D*Xc);
        Xcc = get_block_row(Xc, c, D_range_ext);
        XCp1c = get_block_row(Xc, C+1, D_range_ext);
        Dc = get_block_col(D, c, D_range_ext);
        DCp1 = get_block_col(D, C+1, D_range_ext);

        cost = cost + normF2(Yc - Dc *Xcc - DCp1*XCp1c);
        for i = 1 : C 
            if i == c 
                continue
            end 
            % Di = get_block_col(D, i, D_range_ext);
            Xic = get_block_row(Xc, i, D_range_ext);
            cost = cost + normF2(Xic);
        end 
        cost = 0.5*cost;
    end 
    %%
    function cost = calc_cost(Xc)
       cost = calc_F(Xc) + lambda*norm1(Xc);
    end 
    %%
    DctDc = get_block(DtD, c, c, D_range_ext, D_range_ext);
    DctYc = get_block(DtY, c, c, D_range_ext, Y_range);
    DCp1tDCp1 = get_block(DtD, C+1, C+1, D_range_ext, D_range_ext);
    DCp1tDc = get_block(DtD, C+1, c, D_range_ext, D_range_ext);
    DctDCp1 = DCp1tDc';
    DCp1tYc = get_block(DtY, C+1, c, D_range_ext, Y_range);
    range_c = D_range_ext(c)+1: D_range_ext(c+1);
    range_Cp1 = D_range_ext(C+1) + 1: D_range_ext(C+2);
    DtYc = get_block_col(DtY, c, Y_range);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %================== block:  ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DtYc2 = DtYc;
    DtYc2(range_c,:) = 2*DtYc(range_c,:);
    DtYc2(range_Cp1,:) = 2*DtYc(range_Cp1,:);
    function g = grad(Xc)
        g0 = DtD*Xc - DtYc;
        g1 = Xc;
        Xcc = get_block_row(Xc, c, D_range_ext);
        XCp1c = get_block_row(Xc, C+1, D_range_ext);
        g1(range_c,:) = DctDc*Xcc + DctDCp1*XCp1c - DctYc;
        g1(range_Cp1, :) = DCp1tDCp1*XCp1c + DCp1tDc*Xcc - DCp1tYc;
        g = g0 + g1;% + DtYc2;
    end 
    %%
    function g = num_grad(Xc)
        g = zeros(size(Xc));
        ep = 1e-6;
        for i = 1:size(Xc,1);
            for j = 1: size(Xc,2);
                Xcp = Xc; Xcp(i,j) = Xcp(i,j) + ep;
                Xcm = Xc; Xcm(i,j) = Xcm(i,j) - ep;
                g(i,j) = (calc_F(Xcp) - calc_F(Xcm))/(2*ep);
            end 
        end 
    end 
    %%
    function g = grad2(Xc)
        g0 = DtD*Xc;
%         g1 = Xc;
        Xcc = get_block_row(Xc, c, D_range_ext);
        XCp1c = get_block_row(Xc, C+1, D_range_ext);
        Xc(range_c,:) = DctDc*Xcc + DctDCp1*XCp1c;        
        Xc(range_Cp1, :) = DCp1tDCp1*XCp1c + DCp1tDc*Xcc;
        g = g0 + Xc - DtYc2;
    end
    %% for checking gradient
%     normF2(grad2(Xc) - num_grad(Xc)) 
%     normF2(grad(Xc) - grad(Xc))
%     pause
    %% ========= Main FISTA ==============================
    opts.tol = 1e-8;
    opts.max_iter = 300;
    Xc = fista(@grad2, Xc, L, lambda, opts, @calc_cost);
 
end 