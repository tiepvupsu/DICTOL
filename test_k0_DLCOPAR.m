% test_k0 
dataset = 'myYaleB';
N_train = 10;
k = 8;
lambda1 = 0.001; 
lambda2 = 0.01;
% lambda3 = 0;
max_test = 2;
acc =zeros(9, max_test);
for k0_id = 1:9
	for test_id = 1: max_test
		acc(k0_id, test_id) = ...
			DLCOPAR_top2(dataset, N_train, k, (k0_id -1) *5, lambda1, lambda2);
	end 
end

[(0:10:80)' mean(acc, 2 )]

