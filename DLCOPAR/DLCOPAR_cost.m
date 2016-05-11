function cost = DLCOPAR_cost(Y, Y_range, D, D_range_ext, X, opts)
% function cost = DLCOPAR_cost(Y, Y_range, D, D_range_ext, X, opts)
% Calculating cost of DLCOPAR 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------    
    C = numel(Y_range) - 1;    
    eta    = opts.eta;
    lambda = opts.lambda;
    cost = lambda*norm1(X);
    cost1 = normF2(Y - D*X);
    cost2 = 0;
    DCp1 = get_block_col(D, C+1, D_range_ext);
    DCp1_range = D_range_ext(C+1)+1: D_range_ext(C+2);
    for c = 1: C+1
        Dc = get_block_col(D, c, D_range_ext);
        if c < C+1  
            Dc_range = D_range_ext(c)+1: D_range_ext(c+1);

            Yc = get_block_col(Y, c, Y_range);
            Xc = get_block_col(X, c, Y_range);
            Xcc = Xc(Dc_range,:);
            XCp1c = Xc(DCp1_range,:);
            cost1 = cost1 + normF2(Yc - Dc*Xcc - DCp1*XCp1c);

            Xc(union(Dc_range, DCp1_range),:) = [];
            cost1 = cost1 + normF2(Xc);
        end 
        Dcom_c = D;
        Dcom_c(:, D_range_ext(c)+1: D_range_ext(c+1)) = [];
        cost2 = cost2 + normF2(Dcom_c'*Dc);
    end 
    cost = cost + cost1 + .5*eta*cost2;    
end 