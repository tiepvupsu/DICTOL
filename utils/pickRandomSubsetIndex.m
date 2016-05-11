function [ids, ids_com] = pickRandomSubsetIndex(n, k)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 22:59:40 2016
	% Last modified	: Tue Jan 26 22:59:41 2016
	% Description	: pick a k-element subset of the set 1: n 
	% 	INPUT:
	%		n: number of elements
	%		k: number of picked elements
	% 	OUTPUT: 
	%		ids: vector of picked indices 
	% 		ids_com: vector of unpicked indices
	%
	%% ================== end File info ==========================
    id_mixed = randperm(n);
    ids = id_mixed(1: k);
    if nargout == 2 
    	ids_com = setdiff(1:n, ids);
    end 
end 