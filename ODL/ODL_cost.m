function cost = ODL_cost(Y, D, X, lambda)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 04/07/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	cost = 0.5*normF2(Y - D*X) + lambda*norm1(X);
end 