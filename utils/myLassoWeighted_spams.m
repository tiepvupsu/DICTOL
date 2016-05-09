function X = myLassoWeighted_spams(Y,  D, W, lambdaX)
	param.lambda     = lambdaX; % not more than 20 non-zeros coefficients
	param.lambda2    = 0;
	param.numThreads = 1; % number of processors/cores to use; the default choice is -1
	param.mode       = 2;        % penalized formulation

	Y = normc(Y);
	X = mexLassoWeighted(Y,  D, W, param);
end