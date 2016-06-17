 function cost = D2L2R2_cost(Y, Y_range, D, D_range, X, opts)
    C = numel(Y_range) - 1;
    %% cost of nuclear norm part 
    cost_rank = 0;
    for c = 1: C 
        Dc = get_block_col(D, c, D_range);
        cost_rank = cost_rank + nuclearnorm(Dc);
    end 
    %% total cost
    cost = 0.5*(normF2(Y - D*X) + FDDL_fidelity(Y, Y_range, D, D_range, X)) + ...
          opts.lambda1*norm1(X) + ...
          0.5*opts.lambda2*FDDL_discriminative(X, Y_range) + ...
          opts.alpha*cost_rank;
end 