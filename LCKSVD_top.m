function best_acc = LCKSVD_top(dataset, N_train, k, ...
        sparsitythres, valpha, vbeta)
    % Syntax: LCKSVD_top(dataset, N_train, k, sparsitythres, valpha, vbeta)
    addpath(genpath('utils'));  % add K-SVD box
    addpath(genpath('LCKSVD'));  % add K-SVD box
    addpath('build_spams');
    addpath('utils');
    best_acc = zeros(1, 2);
    if nargin == 0 
        dataset = 'myARgender';
        N_train = 350;
        k = 25;
        dataset = 'myYaleB';
        N_train = 10;
        k = 8;        
        sparsitythres = 10;
        valpha = 0.01;
        vbeta = 0.01;
    end 

    [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
        dataset, N_train);

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
    tic;
    [D1,X1,T1,W1] = labelconsistentksvd1(Y_train,Dinit,Q_train,Tinit,H_train,iterations,sparsitythres,sqrt_alpha);
%     save('.\trainingdata\dictionarydata1.mat','D1','X1','W1','T1');
    rt = toc;
    if strcmp(dataset, 'mySynthetic')
        figure(1);

        display_network(D1);
        drawnow();
    end 

    fprintf('done!');
    %% classification process
    [prediction1,acc] = classification(D1, W1, Y_test, H_test, sparsitythres);
    fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f ', acc);

    fn = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), '_b_', ...
        num2str(vbeta), '_', getTimeStr(), '_1.mat'));
    disp(fn);
   if strcmp(dataset, 'mySynthetic')
        save(fn, 'D1', 'k', 'acc', 'rt');
    else 
        save(fn, 'acc', 'rt');
   end 
    best_acc(1) = acc;

    %% ========= LCKSVD2 ==============================  
    % run LC k-svd training (reconstruction err + class penalty + classifier err)
    fprintf('\nDictionary and classifier learning by LC-KSVD2...');
    tic;
    [D2,X2,T2,W2] = labelconsistentksvd2(Y_train,Dinit,Q_train,Tinit,H_train,Winit,iterations,sparsitythres,sqrt_alpha,sqrt_beta);
    rt = toc;
%     save('.\trainingdata\dictionarydata2.mat','D2','X2','W2','T2');
    fprintf('done!');
    if strcmp(dataset, 'mySynthetic')
        figure(2);

        display_network(D2);
        drawnow();
    end 


    [prediction2,acc] = classification(D2, W2, Y_test, H_test, sparsitythres);
    fprintf('\nFinal recognition rate for LC-KSVD2 is : %.03f\n ', acc);
    fn = fullfile('results', 'LCKSVD', strcat(dataset, '_N_', ...
        num2str(N_train), '_k_', num2str(k), '_a_', num2str(valpha), ...
        '_b_', num2str(vbeta), '_', getTimeStr(), '_2.mat'));
    if strcmp(dataset, 'mySynthetic')
        save(fn, 'D2', 'k', 'acc', 'rt')
    else 
        save(fn, 'acc', 'rt');
    end 
    best_acc(2) = acc;
end 
