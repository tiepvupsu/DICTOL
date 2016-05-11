function [D, X] = FDDL_init(Y, Y_range, opts)
% FDDL Initialization 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/12/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
        addpath(fullfile('..', 'utils'));
        addpath(fullfile('..', 'ODL'));
        addpath(fullfile('..', 'build_spams'));
        d = 10; n = 10; C = 3;
        Y = normc(rand(d, C*n));
        Y_range = n*(0:C);        
        opts.lambda1 = 0.001;
        opts.k = 3;
        opts.max_iter = 20;
    end 
    %%
    C = numel(Y_range) - 1;
    D_range = opts.k*(0:C);
    D = zeros(size(Y,1), opts.k*C);
    for c = 1: C 
        fprintf('%3d ', c);
        if mod(c, 20) == 0
            fprintf('\n');
        end
        range_Dc = get_range(D_range, c);
        Yc = get_block_col(Y, c, Y_range);
        D(:, range_Dc) =    FDDL_INID(Yc , numel(range_Dc), 'pca');
    end 

    X = zeros(size(D, 2), size(Y, 2));
    ini_par.tau         =     opts.lambda1;
    ini_par.lambda      =     opts.lambda2;
    ini_ipts.D          =     D;
    if size(D,1)>size(D,2)
          ini_par.c        =    1.05*eigs(D'*D,1);
    else
          ini_par.c        =    1.05*eigs(D*D',1);
    end
    for ci =  1: C
        fprintf('%3d ', ci);
        if mod(ci, 20) == 0
            fprintf('\n');
        end
        Yc = get_block_col(Y, ci, Y_range);
        ini_ipts.X      =    Yc;
        [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
        X(:, get_range(Y_range, ci)) = ini_opts.A;
    end
        


    % for c = 1: C 
    %     fprintf('Init class %3d\n', c);
    %     range_c = get_range(Y_range, c);
    %     Yc = Y(:, range_c);
    %     range_Dc = get_range(D_range, c);
    %     [Dc, Xcc] = ODL(Yc, opts.k, opts.lambda1, opts, 'spams');
    %     D(:, range_Dc) = Dc;
    %     X(range_Dc, range_c) = Xcc;
    % end 
end 