function cost = COPAR_cost(Y, Y_range, D, D_range_ext, X, opts)
% function cost = COPAR_cost(Y, Y_range, D, D_range_ext, X, opts)
% Calculating cost function of COPAR with parameters lambda and eta are stored in  `opts.lambda` and `opts.rho`.
% `f(D, X) = 0.5*sum_{c=1}^C 05*||Y - DX||_F^2 + 
%               sum_{c=1}^C ( ||Y_c - D_Cp1 X^Cp1_c - D_c X_c^c||F^2 + 
%           sum_{i != c}||X^i_c||_F^2) + lambda*||X||_1 + 
%           0.5*eta*sum_{i \neq c}||Di^T*Dc||_F^2`
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------    
    C = numel(Y_range) - 1;    
    eta    = opts.eta;
    lambda = opts.lambda;
    cost = lambda*norm1(X);
    cost1 = normF2(Y - D*X);
    DCp1 = get_block_col(D, C+1, D_range_ext);
    DCp1_range = D_range_ext(C+1)+1: D_range_ext(C+2);
    for c = 1: C
        Dc = get_block_col(D, c, D_range_ext);
        Dc_range = D_range_ext(c)+1: D_range_ext(c+1);
        Yc = get_block_col(Y, c, Y_range);
        Xc = get_block_col(X, c, Y_range);
        Xcc = Xc(Dc_range,:);
        XCp1c = Xc(DCp1_range,:);
        
        cost1 = cost1 + normF2(Yc - Dc*Xcc - DCp1*XCp1c);

        Xc(union(Dc_range, DCp1_range),:) = [];
        cost1 = cost1 + normF2(Xc);
    end 
    cost = cost + cost1 + .5*eta*DLSI_term(D, D_range_ext);
end 