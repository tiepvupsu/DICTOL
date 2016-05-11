function Mi = get_block_row(M, i, row_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:44:26 AM
	% Last modified	: 1/27/2016 2:44:30 AM
	% Description	: Get a row block of a big matrix X. 
	% 	INPUT
	%		M        : the big matrix. M = [M1 ; M2 ; ... ; MC]
	%		i        : block index 
	%		row_range: a vector store the last index of each block. row_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================
	range = row_range(i) + 1: row_range(i+1);
	Mi = M(range, :);	
end