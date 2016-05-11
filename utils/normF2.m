function res = normF2(A)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 22:33:47 2016
	% Last modified	: Tue Jan 26 22:33:48 2016
	% Description	: square of the Frobenius Norm of A 
	%% ================== end File info ==========================
	res = norm(A, 'fro')^2;
end 