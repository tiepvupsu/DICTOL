function [A, B] = myForm_mult(M, N, P, Q, k)
	% find A, B such that:
	% myForm(M, N, k)*myForm(P, Q, k) = myForm(A, B, k)
	% -----------------------------------------------
	% Author: Tiep Vu, thv102@psu.edu, 4/10/2016
	%         (http://www.personal.psu.edu/thv102/)
	% -----------------------------------------------
	A = M*P;
	B = M*Q + N*P + k*N*Q;
end 