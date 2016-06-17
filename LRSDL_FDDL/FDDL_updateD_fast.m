function [D, it] = FDDL_updateD_fast(Y, Y_range, D, D_range, X, opts)
    
    F = buildMhat(X*X', D_range, D_range);
    E = Y*buildMhat(X, D_range, Y_range)';

    opts.max_iter = 300;
    opts.tol = 1e-8;
    [D, ~] = ODL_updateD(D, E, F, opts);   
end