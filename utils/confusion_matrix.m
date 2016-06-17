function [A, A_pc] = confusion_matrix(pred, gt)
% function A = confusion_matrix(pred, gt)
% calculate confustion matrix from prediction `pred` and ground truth `gt`
% A(i, j) # of truth class i are misclassified as class j 
% A_pc(i, j) percent of truth class i are misclassified as class j 
% example: gt = [1 1 2 2 2 3 3 3 3]
%		pred  = [1 3 2 1 2 1 3 1 3], then 
% A = [ 1 0 1;
% 		1 2 0;
%		1 1 2 ]
% A_pc = [  .5   0  .5;
% 			.33 .66 0;
%			.25 .25 .5 ] sum of each row = 1 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/20/2016 9:59:18 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if nargin == 0
		gt = [1 1 2 2 2 3 3 3 3];
		pred  = [1 3 2 1 2 1 3 2 3];
	end
	C = max(gt);
	A = zeros(C, C);
    %% A 
	for i = 1: numel(gt) 
		A(gt(i), pred(i)) = A(gt(i), pred(i)) + 1;
    end 
    %% A_pc
    s = sum(A, 2);
    A_pc = A./repmat(s, 1, C);
    if nargin == 0
        disp(A) 
        disp(A_pc)
    end 

end 