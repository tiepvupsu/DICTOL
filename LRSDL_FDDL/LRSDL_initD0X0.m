function [D0, X0] = LRSDL_initD0X0(Y, opts, pars)
    % cost = \arg\min .5*\|Y - D0X0\|_F^2 + \lambda1\|X0\|_1 + .5*lambda2*\|X0 - M0\|_F^2 + eta*\|D0\|_F^2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: test module ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SDDL_Dec6_initD0X0()
%     d = 10;
%     N = 100;
%     Y = normc(rand(d, N));
%     opts.lambda1 = 0.001;
%     opts.lambda2 = 0.01;
%     opts.eta = 0.0;
%     pars.max_iter = 20;
%     pars.show = true;       
%     opts.k0 = 10; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: test module ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = opts.lambda3; 
    k0 = opts.k0;
    D0 =  PickDfromY(Y, [0, size(Y,2)], k0);
    X0 = zeros(size(D0, 2), size(Y,2));
    % max_iter = pars.max_iter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: cost function  ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cost = calc_cost(D0, X0)
        cost = .5*normF2(Y - D0*X0) + lambda3*nuclearnorm(D0) + .5*lambda2*normF2(X0 - buildMean(X0)) + lambda1*norm1(X0);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: cost function  ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    iter = 0;
    cost = calc_cost(D0, X0);
    otpsX = opts;
    optsX.max_iter = 20;
    optsD = opts;
    optsD.max_iter = 20;
    while iter < pars.max_iter 
        iter = iter + 1;
        %% ========= update X0 ==============================
        %X0 = \argmin_X0 .5*\|Y - D0X0\|_F^2 + \lambda1\|X0\|_1 + .5*lambda2*\|X0 - M0\|_F^2
        X0 = myLassoWIntrasmall_fista(Y, D0, lambda1, lambda2, X0, optsX);
%         imagesc(X0);
        if pars.verbal
            disp(iter);
            costX0 = calc_cost(D0, X0);
            disp(costX0);
        end 

        %% ========= update D ==============================
        % \arg\min .5*\|Y - D0X0\|_F^2 + eta*\|D0\|_F^2
        D0 = minRankDict0(Y, X0, lambda3, D0, optsD);
        if pars.verbal
            costD0 = calc_cost(D0, X0);
            disp(costD0);
        end
    end
end 