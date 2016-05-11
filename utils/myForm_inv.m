function [P, Q] = myForm_inv(M, N, k)
	% find P, Q such that:
	% myForm(M, N, k)^(-1) = myForm(P, Q, k)
	% -----------------------------------------------
	% Author: Tiep Vu, thv102@psu.edu, 4/10/2016
	%         (http://www.personal.psu.edu/thv102/)
	% -----------------------------------------------
	P = inv(M);
	Q = -inv(M+k*N)*N*P;
end 