 function cost = LRSDL_cost(Y, Y_range, D, D0, D_range, X, X0, opts)
    nClasses = numel(Y_range) - 1;
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = opts.lambda3;        

    if opts.k0 == 0
        Ybar = Y;
    else 
        Ybar = Y - D0*X0;
    end 
    cost = normF2(Ybar - D*X) + FDDL_fidelity(Ybar, Y_range, D, D_range, X);
    cost = .5*cost + lambda1*(norm1([X; X0]));
    cost2 = FDDL_discriminative(X, Y_range) + normF2(X0 - buildMean(X0));
    cost = cost + .5*lambda2*cost2;
    if opts.k0 == 0
        cost = cost + lambda3*nuclearnorm(D0);
    end
end 