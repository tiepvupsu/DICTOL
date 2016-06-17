function Y = shrinkage_rank(D, lambda)
% function Y = shrinkage_rank(D, lambda)
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Wed Jan 27 00:26:45 2016
	% Last modified	: Wed Jan 27 00:26:46 2016
	% Description	: Solving the problem 
	%		Y = argmin_{Y} 0.5*\|D - Y\|_F^2 + lambda\|Y\|_*
	% References:
	% 1. Cai, Jian-Feng, Emmanuel J. Candès, and Zuowei Shen. 
	% 	"A singular value thresholding algorithm for matrix completion." 
	% 	SIAM Journal on Optimization 20.4 (2010): 1956-1982.
	% 	http://arxiv.org/pdf/0810.3286v1.pdf
	% -----------------------------------------------
	% Author: Tiep Vu, thv102@psu.edu, 4/14/2016
	%         (http://www.personal.psu.edu/thv102/)
	% -----------------------------------------------
    if size(D, 1) > size(D, 2)
        [U, S, V] = svd(D, 'econ');
    else 
        [V, S, U] = svd(D', 'econ');
    end
	s = diag(S);
	s1 = max(0, s - lambda);
	S1 = diag(s1);
	Y = U*S1*V';
end