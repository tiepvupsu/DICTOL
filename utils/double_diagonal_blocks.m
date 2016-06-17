function A2 = double_diagonal_blocks(A, row_range, col_range)
% function A = double_diagonal_blocks(A, row_range, col_range)
% Double diagonal blocks of a block matrix A whose row range and col range 
% are `row_range` and `col_range`
% required: `numel(row_range) == numel(col_range)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 % test mode 
		row_range = [0 3 5 8];
		col_range = [0 4 7 10];
		A = rand(row_range(end), col_range(end));
		Ain = A;
		subplot(1, 3, 1); imagesc(A); title('input'); colormap jet;
	end 
	%% check requirements
	if numel(row_range) ~= numel(col_range)
		error('number of blocks in each dimension of the input matrix must be the same');
	end 
	%% MAIN
	mask_diagonal_id = build_diagonal_mask(row_range, col_range);
	A(mask_diagonal_id) = 2*A(mask_diagonal_id);
	%% test mode
	if nargin == 0 
		subplot(1, 3, 2); imagesc(A); title('output');
		subplot(1, 3, 3); imagesc(A - Ain); title('difference');
		A2 = [];
	end 
end 
