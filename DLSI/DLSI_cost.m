function cost = DLSI_cost(Y, Y_range, D, D_range, X, opts)        
    C = numel(Y_range) - 1;
    cost = 0;
    for i = 1: C 
        cost =  opts.lambda*norm1(X{i});
        D_range_i = get_range(D_range, i);        

        Di = D(:, D_range_i);        
        Yi = get_block_col(Y, i, Y_range);
        cost = cost + 0.5*normF2(Yi - Di*X{i}) + norm1(X{i});

        D_com_range_i = setdiff(1:D_range(end), D_range_i);
        D_com_i = D(:, D_com_range_i);
        cost = cost + 0.5*opts.eta*normF2(D_com_i'*Di);
    end 
end 