function acc = calc_acc(pred, ground_truth)
	acc = double(sum(pred == ground_truth))/numel(ground_truth);
end 