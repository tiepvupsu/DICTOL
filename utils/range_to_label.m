function label = range_to_label(range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 3:03:08 AM
	% Last modified	: 1/27/2016 3:03:10 AM
	% Description	: Convert range to label 
	% 		range = [0, 3, 5] then label = [1 1 1 2 2]
	%% ================== end File info ==========================
	C = numel(range) - 1; % number of blocks
	label = ones(1, range(end));
	for i = 2 : C
		n = range(i+1) - range(i);
		label(range(i)+1: range(i + 1)) = i*ones(1, n);
	end
end