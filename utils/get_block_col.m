function Mi = get_block_col(M, i, col_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:39:32 AM
	% Last modified	: 1/27/2016 2:39:36 AM
	% Description	: Get a column block of a big matrix M. 
	% 	INPUT
	%		M        : the big matrix. M = [M1 , M2 , ... , MC]
	%		i        : block index 
	%		col_range: a vector store the last index of each block. col_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================
	range = col_range(i) + 1: col_range(i+1);
	Mi = M(:, range);	
end