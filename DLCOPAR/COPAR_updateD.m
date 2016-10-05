function D = COPAR_updateD(Y, Y_range, D, X, opts) % and DCp1 
% function D = COPAR_updateD(Y, Y_range, D, X, opts) 
% update D in COPAR, including both PARTICULAR dictionaries and the 
% COMMON dictionary.
% The algorithm used here is the efficient algorithm presented in LRSDL paper 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0 
        clc;
        addpath(fullfile('..', 'utils'));
        addpath(fullfile('..', 'DLSI'));
        C = 10;    N = 7;    d = 300;      
        k = 7;
        opts.k0          = 10;
        opts.lambda      = 0.01;
        opts.eta         = 0.1;
        opts.D_range     = k* (0:C);
        opts.D_range_ext = [opts.D_range opts.D_range(end)+opts.k0];
        opts.max_iter    = 10;
        opts.verbose      = true;
        
%         Y       = normc(rand(d, C*N));        
%         Y_range = N * (0:C);
%         D       = normc(rand(d, opts.D_range_ext(end)));
%         X       = 0.01*rand(size(D,2), size(Y,2));
%         save('tmp2.mat', 'Y', 'Y_range', 'D', 'X');

        load('tmp2.mat', 'Y', 'Y_range', 'D', 'X');
    end         
    
    %%
    C              = numel(Y_range) - 1;
    D_range_ext    = opts.D_range_ext;    
    DCp1           = get_block_col(D, C+1, D_range_ext);
    optsD          = opts;
    optsD.verbose   = false;
    optsD.max_iter = 100;    
    Yhat           = zeros(size(Y));
    %% ========= update Dc ==============================
    for c = 1: C 
        %Dc = arg\min_Dc \|Ychat - Dc*Xcc\|_F^2 + \|Ycbar - Dc*Xcc\| + 2*eta\|A*Dc\|_F^2
        % = \arg\min_Dc \| [Ychat Ycbar] - Dc*[Xcc Xcc]\|_F^2 + 2*eta\|A*Dc\|_F^2
        % and solved using DLSI_updateD
        Dc_range = D_range_ext(c)+1: D_range_ext(c+1);
        Yc_range = Y_range(c)+1: Y_range(c+1);
        Yc       = Y(:, Yc_range);
        Dc       = D(:, Dc_range);
        Xc       = X(:, Yc_range);
        Xcc      = get_block_row(Xc, c, D_range_ext);
        XCp1c    = get_block_row(Xc, C+1, D_range_ext);
        Ychat    = Yc - D*Xc + Dc*Xcc;
        Ycbar    = Yc - DCp1*XCp1c;
        E        = (Ychat + Ycbar)*Xcc';
        F        = 2*Xcc*Xcc';
        A        = D;
        A(:,Dc_range)     = []; 
        D(:, Dc_range)    = DLSI_updateD(Dc, E, F, A', opts.eta, optsD);
        Yhat(:, Yc_range) = Yc - D(:, Dc_range)*Xcc;            
    end 
    %% ========= DCp1 ==============================
    XCp1       = X(D_range_ext(C+1) + 1 : D_range_ext(C+2), :);
    Ybar       = Y - D(:, 1: D_range_ext(end-1))*X(1: D_range_ext(end-1), :);
    E          = (Ybar + Yhat)*XCp1';
    F          = 2*XCp1*XCp1';
    A          = D(:, 1: D_range_ext(C+1));
    DCp1_range = D_range_ext(C+1) + 1: D_range_ext(C+2);
    D(:, DCp1_range) = DLSI_updateD(D(:, DCp1_range), E, F, A', opts.eta,  optsD);         
    %%
    if nargin == 0
        D = [];
    end 
end 