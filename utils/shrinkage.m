function X = shrinkage(U, lambda)
	%% =================== Info =========================
	% Author: Tiep Vu
	% Date created: 8/26/2015 7:15:00 PM
	% Description: Solve: X = \arg\min 0.5*||X - U||_F^2 + \lambda ||X||_1
	% Last modified: 8/27/2015 9:46:36 AM
	%% =========================================================
	X = max(0, U - lambda) + min(0, U + lambda);
end
	
