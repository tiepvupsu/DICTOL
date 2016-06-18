function Mi = get_block_row(M, C, row_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:44:26 AM
	% Last modified	: 1/27/2016 2:44:30 AM
	% Description	: Get a row block of a big matrix X. 
	% 	INPUT
	%		M        : the big matrix. M = [M1 ; M2 ; ... ; MC]
	%		C        : block index 
	%		row_range: a vector store the last index of each block. row_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================
	id_sel = [];
	for i = 1: numel(C)
		c = C(i);
		id_sel = [id_sel, row_range(c) + 1: row_range(c+1)];
	end 
	Mi = M(id_sel, : , :);
end