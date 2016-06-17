function res = nuclearnorm(D)
% Return the nuclear norm of the input matrix D    
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/6/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	res = sum(svd(D));
end

% python syntax:
% def nuclearnorm(D):
% 	return sum(np.linalg.svd(D, compute_uv=0))