function test_LRSDL_k0()
	dataset = 'myARgender';
	N_train = 40;
	
	for k0 = 0:10:30 
		for test_id = 1:2

	LRSDL_top(dataset, N_train, k, k0, lambda1, lambda2, lambda3)