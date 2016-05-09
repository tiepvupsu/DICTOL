function W = NNQP_ADMM(F, G, p, b, delta, C, W, opts)
% Solve the problem:
% w = arg min_{w} 0.5*w'*Q*w + p'*w subject to:
% 1. w >= 0
% 2. A*w = b where A = myform(e, 0, C) with e = [1, 1, 1, ...., 1]';
% Q has form myForm(F, G, C); F, G are symmetric nonsingular matrices.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/10/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 % test case 
		addpath(fullfile('..', 'utils'));
		K = 20;
		C = 50;
		F = rand(K, K); F = F'*F;
		G = rand(K, K); G = G'*G;
		p = rand(K*C,1);
		delta = 2;
		W = rand(K, C);
		opts.tol = 1e-8;
		opts.max_iter = 1000;
		%% MATLAB solution for small dimensions.
		Q = myForm(F, G, C);
		e = ones(K,1);
		A = myForm(e', zeros(size(e')), C);
		b = delta*ones(C, 1);
		w = vec(W);
		tic 
		[w_matlab, fval] = quadprog(Q, p, [], [], A, b, zeros(size(w)), [], w);
		toc
		W_matlab = reshape(w_matlab, K, C);
		check_cal = 0; % checking calculation
		profile off;
		profile on;
	end 
	tic
	%%

	function cost = calc_cost(W)
		w = vec(W);
		cost = 0.5*w'*myForm_mult_vec(F, G, w, C) + w'*p;
	end 
	rho = 100;
% 	b = delta*ones(C, 1);
	K = size(G, 1);
	%% Ref: section 7.5 in report
	%% For R 
	% Fhat = R/rho + eye(K); Ghat = G/rho.
	if check_cal
		QQ = [Q/rho + eye(K*C), A'/rho; A zeros(size(A,1), size(A,1))];
	end 
	[Fbar, Gbar] = myForm_inv(F/rho + eye(K), G/rho, C);
	if check_cal
		fprintf('check R\n');
		R_def = inv(Q/rho + eye(K*C));
		R_cal = myForm(Fbar, Gbar, C);
		norm(R_def - R_cal)
	end 
	% pause
	%% For S 
	% e'*F*e = sum(vec(F));
	fbar = sum(Fbar,2);
	gbar = sum(Gbar,2);
	[S1, S2] = myForm_inv(sum(fbar), sum(gbar), C);
	if check_cal
		fprintf('check S\n');
		S_def = inv(A*R_def*A');
		S_cal = myForm(S1, S2, C);
		norm(S_def - S_cal)
	end 
	%% RAtS 
	[RAtS1, RAtS2] = myForm_mult(fbar, gbar, S1, S2, C);
	%% RAtSAR
	[RAtSAR1, RatSAR2] = myForm_mult(RAtS1, RAtS2, fbar', gbar', C);
	N1 = Fbar - RAtSAR1;
	N2 = Gbar - RatSAR2;
	%%
	if check_cal
		N = myForm(N1, N2, C);
		RAtS = myForm(RAtS1, RAtS2, C);
		R = myForm(Fbar, Gbar, C);
		S = myForm(S1, S2, C);
		invQQ2 = [N RAtS; rho*S*A*R -rho*S];
	end 
	
	% (inv(QQ) - invQQ2)
	% inv(Q/rho)
	% myForm(Fbar, Gbar, C) - inv(Q/rho + eye(K*C))
	% pause

	iter = 0;
	U = zeros(size(W));
	Z = U;
	P = reshape(p, K, C);
	if nargin == 0 % check inversion
		QQ = [Q/rho + eye(K*C), A'/rho; A zeros(size(A,1), size(A,1))];
		QQ1 = inv(QQ);
		myForm(N1, N2, C);
% 		pause 
	end 
	%%
	while iter < opts.max_iter 
		iter = iter + 1;
        
		%% update W 
		V = Z - U - P/rho;
		w = myForm_mult_vec(N1, N2, vec(V), C) + myForm_mult_vec(RAtS1, RAtS2, b, C);
		
% 		w1 = inv(QQ)*[vec(V); b];
% 		w = w1(1: K*C);
        W = reshape(w, K, C);
		%% update Z 
		Z_new = max(0, W + U);
        rk = norm(W - Z);
        sk = rho*norm(Z_new - Z);
        if rk < opts.tol && sk < opts.tol
            break;
        end 
        Z = Z_new;
		%% update U 
		U = U + W - Z;
		%%
	end 
	toc
	if nargin == 0
        sum(W,1)
        sum(W_matlab,1)
%         [vec(W) vec(W_matlab)]
		norm(W - W_matlab)
		calc_cost(W) - calc_cost(W_matlab)
% 		profile viewer;
		pause;
	end
end 
