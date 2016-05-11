function y = myForm_mult_vec(M, N, x, k)
	% y = myForm(M, N, k)*vec(X)
	% -----------------------------------------------
	% Author: Tiep Vu, thv102@psu.edu, 4/10/2016
	%         (http://www.personal.psu.edu/thv102/)
	% -----------------------------------------------
	X = reshape(x, size(M,2), k);
	y = vec(M*X) + repmat(sum(N*X, 2), k, 1);
end 