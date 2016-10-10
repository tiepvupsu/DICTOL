function [acc, rt] = LCKSVD_wrapper(Y_train, label_train, Y_test, label_test,...
                    k, sparsitythres, valpha, vbeta)
% function [acc, rt] = LCKSVD_wrapper(Y_train, label_train, Y_test, label_test,...
                    % k, sparsitythres, valpha, vbeta)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/19/2016 2:47:38 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
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
        [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
            dataset, N_train);
    end 
    acc = zeros(1, 2);
    rt = zeros(1, 2);
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
    [Dinit,Tinit,Winit,Q_train] = initialization4LCKSVD(Y_train, H_train, ...
        dictsize, iterations4ini, sparsitythres);
    fprintf('done!');
    %% ========= LCKSVD1 ==============================  
    % run LC K-SVD Training (reconstruction err + class penalty)
    fprintf('\nDictionary and classifier learning by LC-KSVD1...');
    tic;
    [D1,X1,T1,W1] = labelconsistentksvd1(Y_train, Dinit, Q_train, Tinit, ...
        H_train,iterations,sparsitythres,sqrt_alpha);
    rt(1) = toc;
    fprintf('done!');
    %% classification process
    [prediction1, acc(1)] = classification(D1, W1, Y_test, H_test, sparsitythres);
%     fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f \n', acc);
    %% ========= LCKSVD2 ==============================  
    % run LC k-svd training (reconstruction err + class penalty + classifier err)
    fprintf('\nDictionary and classifier learning by LC-KSVD2...');
    tic;
    [D2,X2,T2,W2] = labelconsistentksvd2(Y_train, Dinit, Q_train, ...
        Tinit, H_train, Winit, iterations, sparsitythres, sqrt_alpha, sqrt_beta);
    rt(2) = toc;
    fprintf('done!\n');
    [prediction2, acc(2)] = classification(D2, W2, Y_test, H_test, sparsitythres);
end 