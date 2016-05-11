function Mij = get_block(M, i, j, row_range, col_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:46:04 AM
	% Last modified	: 1/27/2016 2:46:09 AM
	% Description	: Get a block matrix of a big matrix X. 
	% 	INPUT
	%		M        : the big matrix. 
	%			M = [	M11, M12, ..., M1m]
	% 					M21, M22, ..., M2m]
	%					..................
	% 					Mn1, Mn2, ..., Mnm]
	%		i        : row block index 
	%		j        : column block index 
	%		row_range: a vector storing the last index of each block. row_range(1) = 0.
	%					i-th block is indexed by row_range(i)+1: row_range(i+1).
	%		col_range: a vector storing the last index of each block. row_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================
	range1 = row_range(i) + 1: row_range(i + 1);
	range2 = col_range(j) + 1: col_range(j + 1);
	Mij = M(range1, range2);
end