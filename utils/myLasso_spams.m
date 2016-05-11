function X = lasso_spams(Y, D, lambdaX)
	param.lambda     = lambdaX; % not more than 20 non-zeros coefficients
	param.lambda2    = 0;
	param.numThreads = 1; % number of processors/cores to use; the default choice is -1
	param.mode       = 2;        % penalized formulation


	%% ========= normalize Y ==============================
	Y = normc(Y);
	X = mexLasso(Y, D, param);

	
end