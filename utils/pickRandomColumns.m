function [Yk, Yk_com] = pickRandomColumns(Y, k)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 23:11:01 2016
	% Last modified	: Tue Jan 26 23:11:02 2016
	% Description	: Pick randomly k columns of Y, k must be less than size(Y,2)
	% 	INPUT:
	%		Y: a matrix 
	%		k: number of picked columns 
	% 	OUTPUT: 
	%		Yk: a submatrix comprising the picked columns of Y 
	%		Yk_com: a complementary matrix (unpicked columns)
	%% ================== end File info ==========================
	N = size(Y,2)
	if k > N 
		error('k must be less than #columns of Y\n');
	end 
	[ids, ids_com] = pickRandomSubsetIndex(N, k);
	Yk = Y(:, ids);
	if nargout == 2 
		Yk_com = Y(:, ids_com);
	end 
end 