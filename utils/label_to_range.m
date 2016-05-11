function range = label_to_range(label)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 3:00:12 AM
	% Last modified	: 1/27/2016 3:00:16 AM
	% Description	: Convert label to range, label = [1 1 1 2 2 2 2 3 3]
	%                 then range = [0, 3, 7, 9]
	%		
	%
	%% ================== end File info ==========================
	C = label(end); % number of blocks
	range = zeros(1, max(label) +1);
	t = 0;
	for i = 1 : C
		t = t + numel(find(label == i));
		range(i + 1) = t;
	end
end 


