function [Y_trn, label_trn, Y_tst, label_tst] = picktrntst(Y, Y_range, N_trn_c)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: Tue Jan 26 23:43:24 2016
	% Last modified	: Wed Jan 27 00:08:00 2016
	% Description	: pick training set and test set from data
	% 	INPUT:
	%		Y: data (each column is an observation)
	%		label: label of data (start from 1 to C - nmber of class)
	%		N_trn_c: number of training samples in each class
	% 	OUTPUT: 
	%		Y_trn: picked training data 
	%		label_trn: label of training data 
	%		Y_tst: test data 
	%		label_tst: test label
	%% ================== end File info ==========================

	%% ========= Main code ==============================
	% [Y, label] = normalizeDataLabel(Y, label);

	% Y_range = label_to_range(label);
	C       = numel(Y_range) - 1; % number of classes 
    d = size(Y,1);
	N_all = Y_range(C+1);
	N_trn = C*N_trn_c;
	N_tst = N_all - N_trn;

	Y_trn     = zeros(d, N_trn);
	Y_tst     = zeros(d, N_tst);
	label_trn = zeros(1, N_trn);
	label_tst = zeros(1, N_tst);

	cur_trn = 0;
	cur_tst = 0;
    for c = 1: C 
        Yc = get_block_col(Y, c, Y_range);

		N_all_c = size(Yc, 2);
		N_tst_c = N_all_c - N_trn_c;
		label_trn(:, cur_trn + 1: cur_trn + N_trn_c) = c*ones(1, N_trn_c);
		label_tst(:, cur_tst + 1: cur_tst + N_tst_c) = c*ones(1, N_tst_c);

        idx = randperm(N_all_c);
		
		Y_trn(:, cur_trn + 1: cur_trn + N_trn_c) = Yc(:, idx(1: N_trn_c));
		Y_tst(:, cur_tst + 1: cur_tst + N_tst_c) = Yc(:, idx(N_trn_c + 1: end));

		cur_trn = cur_trn + N_trn_c;
		cur_tst = cur_tst + N_tst_c;
    end 

end 

