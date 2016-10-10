function Xc = COPAR_updateXc(DtD, DtY,  Y_range, Xc, c, L, opts) 
% function Xc = COPAR_updateXc(DtD, DtY,  Y_range, Xc, c, L, opts) 
% * Update Xc in COPAR (page 189-190 COPAR)
% see COPAR paper: 
% http://www.cs.zju.edu.cn/people/wangdh/papers/draft_ECCV12_particularity.pdf
%  cost = normF2(Yc - D*Xc) + normF2(Yc - DcXcc - DCp1*XCp1c) + 
%           sum_{i \neq c, 1 \leq i \leq C} normF2(Xic);
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/12/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
        addpath(fullfile('..', 'utils'));
        fprintf('Test most\n');
        C = 3;    N = 10;    d = 30;
        k = 10;
        k0 = 3;    
        c = 2;
%         Y       = normc(rand(d, C*N));
%         D       = normc(rand(d, C*k + k0));
%         DtD     = D'*D;
%         DtY     = D'*Y;
%         Y_range = N*(0:C);
%         Yc      = get_block_col(Y, c, Y_range);
%         Xc      = zeros(size(D,2), size(Yc,2));
%         L       = 2*max(eig(DtD))+10;
%         D_range = k* (0:C);
%         D_range_ext = [D_range D_range(end)+k0];
%         save('tmp.mat', 'Y', 'D', 'Y_range', 'D_range_ext', 'L');
        load('tmp.mat', 'Y', 'D', 'Y_range', 'D_range_ext', 'L');
        
        DtD = D'*D;
        DtY = D'*Y;
        Yc      = get_block_col(Y, c, Y_range);
        Xc      = zeros(size(D,2), size(Yc,2));
        D_range = D_range_ext(1: end-1);
        
        opts.D_range     = k* (0:C);
        opts.D_range_ext = D_range_ext;
        opts.k0          = k0;
        opts.lambda      = 0.01;
        opts.eta         = 0.1;
        opts.max_iter    = 300;
        opts.verbose      = true;
        opts.check_grad  = true;  
        
    end 
    %%
    C           = numel(opts.D_range_ext) - 2;
    D_range_ext = opts.D_range_ext;  
    lambda      = opts.lambda/2;
    %%
    function cost = calc_f(Xc)
        cost  = normF2(Yc - D*Xc);
        Xcc   = get_block_row(Xc, c, D_range_ext);
        XCp1c = get_block_row(Xc, C+1, D_range_ext);
        Dc    = get_block_col(D, c, D_range_ext);
        DCp1  = get_block_col(D, C+1, D_range_ext);
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
    function cost = calc_F(Xc)
       cost = calc_f(Xc) + lambda*norm1(Xc);
    end 
    %%
    DctDc              = get_block(DtD, c, c, D_range_ext, D_range_ext);
    DCp1tDCp1          = get_block(DtD, C+1, C+1, D_range_ext, D_range_ext);
    DCp1tDc            = get_block(DtD, C+1, c, D_range_ext, D_range_ext);
    DctDCp1            = DCp1tDc';
    range_c            = D_range_ext(c)+1: D_range_ext(c+1);
    range_Cp1          = D_range_ext(C+1) + 1: D_range_ext(C+2);
    DtYc               = get_block_col(DtY, c, Y_range);
    DtYc2              = DtYc;
    DtYc2(range_c,:)   = 2*DtYc(range_c,:);
    DtYc2(range_Cp1,:) = 2*DtYc(range_Cp1,:);
    %%
    function g = grad(Xc)
        g0 = DtD*Xc;
        Xcc              = get_block_row(Xc, c, D_range_ext);
        XCp1c            = get_block_row(Xc, C+1, D_range_ext);
        Xc(range_c,:)    = DctDc*Xcc + DctDCp1*XCp1c;        
        Xc(range_Cp1, :) = DCp1tDCp1*XCp1c + DCp1tDc*Xcc;
        g = g0 + Xc - DtYc2;
    end
    %% check grad 
    if opts.check_grad 
        check_grad(@calc_f, @grad, rand(size(Xc)));
    end 
    %% ========= Main FISTA ==============================
    L2 = max(eig(DtD)) + max(eig(DctDc')) + max(eig(DCp1tDCp1'));

    opts.tol = 1e-8;
    opts.max_iter = 300;
    Xc1 = fista(@grad, Xc, L, lambda, opts, @calc_F);
    pause
    Xc2 = fista(@grad, Xc, L2, lambda, opts, @calc_F);
    %%
    if nargin == 0 
        Xc = [];
    end 
end 