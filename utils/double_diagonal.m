function A = double_diagonal(A)
% function A = double_diagonal(A)
% dDouble diagonal of a matrix A
% required: `size(A, 1) == size(A, 2)`
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0 % test mode 
		A = rand(4, 4);
		Ain = A;
		subplot(1, 3, 1); imagesc(A); title('input'); colormap jet;
	end 
	%% check requirements
	if size(A, 1) ~= size(A, 2)
		error('The input matrix is not square!!');
	end 
	%% MAIN
	n = size(A, 1);
    A(1:n+1:n^2) = 2*diag(A);
    %%
    if nargin == 0 % test mode
    	subplot(1, 3, 2); imagesc(A); title('output');
    	subplot(1, 3, 3); imagesc(Ain - A); title('difference');
    	A = [];
    end 
end 
