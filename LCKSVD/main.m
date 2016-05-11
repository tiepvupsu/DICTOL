% =========================================================================
% An example code for the algorithm proposed in
%
%   Zhuolin Jiang, Zhe Lin, Larry S. Davis.
%   "Learning A Discriminative Dictionary for Sparse Coding via Label 
%    Consistent K-SVD", CVPR 2011.
%
% Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
% Date: 10-16-2011
% =========================================================================


clear all;
clc;
addpath(genpath('.\ksvdbox'));  % add K-SVD box
addpath(genpath('.\OMPbox')); % add sparse coding algorithem OMP
load('.\trainingdata\featurevectors.mat','training_feats', 'testing_feats', 'H_train', 'H_test');

%% constant
sparsitythres = 30; % sparsity prior
sqrt_alpha = 4; % weights for label constraint term
sqrt_beta = 2; % weights for classification err term
dictsize = 570; % dictionary size
iterations = 50; % iteration number
iterations4ini = 20; % iteration number for initialization

%% dictionary learning process
% get initial dictionary Dinit and Winit
fprintf('\nLC-KSVD initialization... ');
[Dinit,Tinit,Winit,Q_train] = initialization4LCKSVD(training_feats,H_train,dictsize,iterations4ini,sparsitythres);
fprintf('done!');

% run LC K-SVD Training (reconstruction err + class penalty)
fprintf('\nDictionary learning by LC-KSVD1...');
[D1,X1,T1,W1] = labelconsistentksvd1(training_feats,Dinit,Q_train,Tinit,H_train,iterations,sparsitythres,sqrt_alpha);
save('.\trainingdata\dictionarydata1.mat','D1','X1','W1','T1');
fprintf('done!');

% run LC k-svd training (reconstruction err + class penalty + classifier err)
fprintf('\nDictionary and classifier learning by LC-KSVD2...')
[D2,X2,T2,W2] = labelconsistentksvd2(training_feats,Dinit,Q_train,Tinit,H_train,Winit,iterations,sparsitythres,sqrt_alpha,sqrt_beta);
save('.\trainingdata\dictionarydata2.mat','D2','X2','W2','T2');
fprintf('done!');

%% classification process
[prediction1,accuracy1] = classification(D1, W1, testing_feats, H_test, sparsitythres);
fprintf('\nFinal recognition rate for LC-KSVD1 is : %.03f ', accuracy1);

[prediction2,accuracy2] = classification(D2, W2, testing_feats, H_test, sparsitythres);
fprintf('\nFinal recognition rate for LC-KSVD2 is : %.03f ', accuracy2);