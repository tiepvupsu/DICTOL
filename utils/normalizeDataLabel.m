function [normalized_Y, normalized_label] = normalizeDataLabel(Y, label)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 23:46:51 2016
	% Last modified	: Tue Jan 26 23:46:52 2016
	% Description	: normalize each column of Y to norm = 1, label is in form: 
	% 		1,1,..., 1, 2, 2,..., 2, 3, 3..., 3, ... (same-class samples are adjacent) (sorted)
	%% ================== end File info ==========================	
	Y = normc(Y);
	if issorted(label) % if label is sorted 
		normalized_Y = Y;
		normalized_label = label;
	else 
		[normalized_label, p] = sort(label);
		normalized_Y = Y(:, p);
	end 
end 