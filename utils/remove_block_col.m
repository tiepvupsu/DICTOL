function M = remove_block_col(M, c, M_range)
	range = M_range(c)+1: M_range(c+1);
	M(:, range) = [];
end