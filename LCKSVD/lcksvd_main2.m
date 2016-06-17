function LCKVSD_top(dataset, N_train, k, sparsitythres, valpha, vbeta)
    myinit();
    [Y_train, label_train, Y_test, label_test] = pickTrainTest_2(dataset, N_train);



    %%
    C = max(label_train);
    
    H_train = lcksvd_buildH(label_train);
    H_test = lcksvd_buildH(label_test);


    % sparsitythres = 30; % sparsity prior
    sqrt_alpha = sqrt(valpha); % weights for label constraint term
    sqrt_beta = sqrt(vbeta); % weights for classification err term
    dictsize = C*k; % dictionary size
    iterations = 50; % iteration number
    iterations4ini = 20; % iteration number for initialization

    %% dictionary learning process
    % get initial dictionary Dinit and Winit
    fprintf('\nLC-KSVD initialization... ');
    [Dinit,Tinit,Winit,Q_train] = initialization4LCKSVD(Y_train,H_train,dictsize,iterations4ini,sparsitythres);
    fprintf('done!');

    %% ========= LCKSVD1 ==============================  
    % run LC K-SVD Training (reconstruction err + class penalty)
    fprintf('\nDictionary learning by LC-KSVD1...');
    [D1,X1,T1,W1] = labelconsistentksvd1(Y_train,Dinit,Q_train,Tinit,H_train,iterations,sparsitythres,sqrt_alpha);
    save('.\trainingdata\dictionarydata1.mat','D1','X1','W1','T1');
    fprintf('done!');
    %% classification process
    [prediction1,acc] = classification(D1, W1, Y_test, H_test, sparsitythres);
    fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f ', accuracy1);

    fn = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), '_b_', num2str(vbeta), '_1.mat'));
    save(fn, 'acc');


    %% ========= LCKSVD2 ==============================  
    % run LC k-svd training (reconstruction err + class penalty + classifier err)
    fprintf('\nDictionary and classifier learning by LC-KSVD2...')
    [D2,X2,T2,W2] = labelconsistentksvd2(Y_train,Dinit,Q_train,Tinit,H_train,Winit,iterations,sparsitythres,sqrt_alpha,sqrt_beta);
    save('.\trainingdata\dictionarydata2.mat','D2','X2','W2','T2');
    fprintf('done!');


    [prediction2,acc] = classification(D2, W2, Y_test, H_test, sparsitythres);
    fprintf('\nFinal recognition rate for LC-KSVD2 is : %.03f ', accuracy2);
    fn = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), '_b_', num2str(vbeta), '_2.mat'));
    save(fn, 'acc');
end 