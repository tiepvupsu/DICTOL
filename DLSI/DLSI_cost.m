function cost = DLSI_cost(Y, Y_range, D, D_range, X, opts)        
% function cost = DLSI_cost(Y, Y_range, D, D_range, X, opts)        
% Calculating cost function of DLSI with parameters lambda and eta are stored in 
% `opts.lambda` and `opts.rho`
% f(D, X) = 0.5*sum_{c=1}^C 05*||Yc - Dc*Xc||_F^2 + lambda*||X||_1 + 
%           05*eta*sum_{i \neq c} ||Di^T*Dc||_F^2
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    C = numel(Y_range) - 1;
    cost = 0;
    for i = 1: C 
        cost = cost + opts.lambda*norm1(X{i});
        D_range_i = get_range(D_range, i);        
        Di = D(:, D_range_i);        
        Yi = get_block_col(Y, i, Y_range);
        cost = cost + 0.5*normF2(Yi - Di*X{i});
    end 
    cost = cost + .5*opts.eta*DLSI_term(D, D_range);
end 