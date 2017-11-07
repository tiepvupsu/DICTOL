function D = min_rank_dict(Dinit, E, F, lambdaD, opts)
% function D = min_rank_dict(D, E, F, lambdaD, opts)
% This function solves the following problem:
% [D, X] = argmin_D -2trace(ED') + trace(FD'D) + lambdaD ||D||_*
% s.t. ||d_i||_2^2 <= 1, for all i 
% using ADMM:
% INPUT: 
% 	Y: Data 
% 	Dinit: intitial D 
% 	X: sparse code  
% 	lambdaD: regularization term 
% OUTPUT: 
% 	D: 
% -------------- Detail ADMM procedure ---------
% Choose a rho.
% Algorithm summary
% ADMM: D,J = argmin_{D, J} -2trace(ED') + trace(FD'D) + lambdaD||J||_*
% s.t ||d_i||_2^2 <= 1 and J = D
% Alternatively solving:
% (1): D^{k+1} = argmin_D -2trace(ED') + trace(FD'D) + rho/2 ||J - D + U^k||_F^2 
% 		s.t. ||d_i||_2^2 <= 1
% 	this problem can be soved by ODL
% (2): J^{k+1} = argminJ lambdaD||J||_* + rho/2||J - D^{k+1} + U^k||
% 	Solution: shrinkage_rank(D^{k+1} - U^k, lambdaD/rho)
% (3): Update U: U^{k+1} = U^k + J^{k+1} - D^{k+1}
% Stoping cretia:
% ||r^k||_F^2 <= tol, ||s^k||_F^2 <= tol 
% r^k = J^k - D^k 
% s^k = rho(J^{k+1} - J^k) 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/12/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 
        clc;
        addpath(fullfile('..', 'ODL'));
		d = 504;
		k = 10;
		Dinit = normc(rand(d, k));
		N = 380;
		Y = normc(rand(d, N));
		X = rand(k, N);
		E = Y*X';
		F = X*X';
		lambdaD = 0.001;
		opts.verbose = 0;
		opts.max_iter = 1000;
    end 
    %%
   
    if lambdaD == 0
    	D = ODL_updateD(Dinit, E, F, opts);
        if nargin == 0
            D = [];
        end 
    	return;        
    end 
    %%
    function cost = calc_cost(D)
        cost = -2*trace(E*D') + trace(F*D'*D) + lambdaD*nuclearnorm(D);
    end 
    %%
    opts = initOpts(opts);
    rho = 0.25;
    D_old = Dinit ;
    J_old = Dinit;
    U_old = zeros(size(Dinit));
    k = 0;
    I = eye(size(Dinit,2));
    tau = 2.0;
    mu = 10.0;
    optsD = opts;
    optsD.max_iter = 50;
    opts.tol = 1e-5;
    t1 = 0;
    t2 = 0;
    while k < opts.max_iter
	k = k + 1;		
	%% ========= update D ==============================
	% D = argmin_D -2trace(ED') + trace(FD'D) + rho/2 ||J - D + U||_F^2 
	% s.t. ||d_i||_2^2 <= 1
	E1 = E + rho/2*(J_old + U_old);
	F1 = F + rho/2*I;
	D_new = ODL_updateD(D_old, E1, F1, optsD);
	%% ========= update J ==============================
	J_new = real(shrinkage_rank(D_new - U_old, lambdaD/rho));
	%% ========= update U ==============================
	U_new = U_old + J_new - D_new;        
	%% ========= check stop ==============================
	r = J_new - D_new;
	s = rho*(J_new - J_old);
	r_eps = norm(r, 'fro');
	s_eps = norm(s, 'fro')		;
	if (r_eps < opts.tol && s_eps < opts.tol)
		break;
        end
        if opts.verbose
	    cost = calc_cost(D_new);
            fprintf('iter = %2d | cost = %5.4f | r_eps %4.4f | s_eps = %4.4f\n',...
                   k, cost, r_eps, s_eps);
        end 
        D_old = D_new;
        J_old = J_new;
        U_old = U_new;
        if r_eps > mu*s_eps
            rho = rho*tau;
        elseif s_eps > mu*r_eps
            rho = rho/tau;
        end
    end 
    %% ========= return ==============================
    D = D_new;
    if nargin == 0
        D = [];
    end
end
			
		
