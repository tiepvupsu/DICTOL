function [Y_trn, label_trn, Y_tst, label_tst] = picktrntst_wrapper(dataset, N_trn)
	data_fn = fullfile('data', strcat(dataset, '.mat'));
    load(data_fn);
    if ~exist('Y_range', 'var')
        Y_range = label_to_range(label);
    end

    Y = normc(Y);

    [Y_trn, label_trn, Y_tst, label_tst] = picktrntst(Y, Y_range, N_trn);
end
