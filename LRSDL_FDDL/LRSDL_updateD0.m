% function D0 = LRSDL_updateD0(Y, D, D0, Y_range, D_range, opts, X, X0, pars)
function D0 = LRSDL_updateD0(Y, Y_range, D, D_range, D0, X, X0, opts)
    % D0 = armin_D0 0.5*normF2(L - D0X0) + eta/2 * nuclearnorm(D0)
%     Ybar = Y - D*X;
% 
%     mask = build_mask(D_range, Y_range);
%     Yhat = Y - D*(X.*mask);
% 
%     % Yhat = LRSDL_buildYhat(Y, D, X, D_range, Y_range);
% 
%     L = (Ybar + Yhat)/2;
    % D0 = minRankDict(Ysum*X0', D0, X0*X0', eta);    
    L = Y - D*buildMhat(X, D_range, Y_range)/2;
%     norm(L - L2)
%     D0 = minRankDict0(L, X0, opts.lambda3/2, D0, pars);
    D0 = min_rank_dict(D0, L*X0', X0*X0', opts.lambda3, opts);
end