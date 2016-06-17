function pred = DLSI_pred(Y, D, opts)
% function pred = DLSI_pred(Y, D, opts)
% j = \arg\min_j R(y, Dj) with R(y,D) = 0.5*\|y - Dx\|_2^2 + lambda*\|x\|_1;
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    D_range = opts.D_range;
    C = numel(opts.D_range)-1;
    E = zeros(C, size(Y,2));
    opts.show = false;
    opts.max_iter = 300;
    for c = 1: C 
        Dc = get_block_col(D, c, D_range);
        Xc = lasso_fista(Y, Dc, [], opts.lambda, opts);
        R1 = Y - Dc*Xc;
        E(c,:) = 0.5* sum(R1.^2, 1) + opts.lambda*sum(abs(Xc),1);
    end 
    [~, pred] = min(E);
end 
