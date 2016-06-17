% function D0 = LRSDL_updateD0(Y, D, D0, Y_range, D_range, opts, X, X0, pars)
function D0 = LRSDL_updateD0(Y, Y_range, D, D_range, D0, X, X0, opts)

    L = Y - D*buildMhat(X, D_range, Y_range)/2;
    D0 = min_rank_dict(D0, L*X0', X0*X0', opts.lambda3, opts);
end