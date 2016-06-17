function X = myOMP(Y, D, L)
	%% =================== Info =========================
	% Author:  Tiep Vu thv102@psu.edu
	% Date created: 8/27/2015 9:09:35 AM
	% Description: Solve the l0 problem using OMP. This is actuaaly 
	% a wrapper of mexOMP in SPAMS toolbox
	%% =========================================================
	% param.L          = L; % not more than 10 non-zeros coefficients
	% param.eps        = 0.001; % squared norm of the residual should be less than 0.1
	% param.numThreads = 1; % number of processors/cores to use; the default choice is -1

	% Y = normc(Y);
	% X = mexOMP(Y, D, param);
	X = omp(D'*Y, D'*D, L);
end
