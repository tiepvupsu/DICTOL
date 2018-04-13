clc;
clear all;
close all;

addpath(genpath('utils'));  % add K-SVD box
addpath(genpath('LCKSVD'));  % add K-SVD box
addpath('COPAR');
addpath('DLSI');
addpath('LRSDL_FDDL');
addpath('build_spams');
addpath('utils');
addpath('SRC');
addpath('ODL');

dataset = 'myYaleB';
% dataset = 'myAR';
N_train = 10;        
fprintf('Start time\n');
t = getTimeStr();
[dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
    dataset, N_train);

range_train = label_to_range(label_train);
range_test = label_to_range(label_test);


fprintf('================= SRC ====================\n');
opts.lambda = 0.001;
lambda = 0.001;
% acc_src = SRC_top;
acc_src = SRC_wrapper(Y_train, range_train, Y_test, range_test, lambda);
disp(acc_src);
% fprintf('================= LCKSVD =================\n');
% 
% valpha = 0.002;
% vbeta = 0.004;
% k = 10;        
% sparsitythres = 10;
% [acc_lcksvd, rt] = LCKSVD_wrapper(Y_train, label_train, Y_test, label_test,...
%                     k, sparsitythres, valpha, vbeta);
% disp(acc_lcksvd);
% fprintf('\n================= DLSI ===================\n');
% k = 10;
% lambda = 0.001;
% eta = 0.01;
% [acc_dlsi, rt] = DLSI_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, lambda, eta);
% disp(acc_dlsi);
fprintf('================= FDDL ===================\n');
k = 10;
lambda1 = 0.001;
lambda2 = 0.05;
[acc_fddl, rt] = FDDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, lambda1, lambda2);
                        
disp(acc_fddl);
acc_fddl = FDDL_top;
% fprintf('================= COPAR ================\n');
% % acc_COPAR = COPAR_top;
% k=10;
% k0 = 5;
% lambda = 0.001;
% eta = 0.01;
% [acc_COPAR, rt] = COPAR_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda, eta);
% disp(acc_COPAR);
% fprintf('================= D2L2R2 =================\n');
% acc_d2l2r2 = D2L2R2_top;
fprintf('================= LRSDL ==================\n');
% acc_lrsdl = LRSDL_top;
k = 10;
k0 = 5;
lambda1 = 0.001;
lambda2 = 0.01;
lambda3 = 0.02;
[acc_lrsdl, rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda1, lambda2, lambda3);
disp(acc_lrsdl);
fprintf('================= Summaray =================\n');
fprintf('+--------------------------+\n')
fprintf('|  Method    |   Accuray   |\n')
fprintf('+------------+-------------+\n')
fprintf('|  SRC       |   %2.2f%%    |\n', 100*acc_src);
% fprintf('|  LCKSVD1   |   %2.2f%%    |\n', 100*acc_lcksvd(1));
% fprintf('|  LCKSVD2   |   %2.2f%%    |\n', 100*acc_lcksvd(2));
% fprintf('|  DLSI      |   %2.2f%%    |\n', 100*acc_dlsi);
fprintf('|  FDDL      |   %2.2f%%    |\n', 100*acc_fddl);
% fprintf('|  COPAR     |   %2.2f%%    |\n', 100*acc_COPAR);
% fprintf('|  D2L2R2    |   %2.2f%%    |\n', 100*acc_d2l2r2);
fprintf('|  LRSDL     |   %2.2f%%    |\n', 100*acc_lrsdl);
fprintf('+--------------------------+\n')
t = getTimeStr();
fprintf('Finish at \n');
