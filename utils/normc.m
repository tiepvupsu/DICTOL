function res = normc(A)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 22:15:40 2016
	% Last modified	: Tue Jan 26 22:15:42 2016
	% Description	: normalize columns of a matrix (each column has Euclidean norm = 1)
	%		This is a built-in function in some new version of MATLAB
	%% ================== end File info ==========================
	B = A.^2;
	C = sqrt(sum(B));
	res = A./repmat(C, size(A,1), 1);
end