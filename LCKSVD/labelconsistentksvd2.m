% ========================================================================
% Label consistent KSVD algorithm 2
% USAGE: [D,X,T,W]=labelconsistentksvd2(Y,Dinit,Q_train,Tinit,H_train,....
%                      Winit,iterations,sparsitythres,sqrt_alpha,sqrt_beta)
% Inputs
%       Y               -training features
%       Dinit           -initialized dictionary
%       Q_train         -optimal code matrix for training feature 
%       Tinit           -initialized transform matrix
%       H_train         -labels matrix for training feature 
%       Winit           -initialized classifier parameters
%       iterations      -iterations for KSVD
%       sparsitythres   -sparsity threshold for KSVD
%       sqrt_alpha      -contribution factor
%       sqrt_beta       -contribution factor
% Outputs
%       D               -learned dictionary
%       X               -sparsed codes
%       T               -learned transform matrix
%       W               -learned classifier parameters
%
% Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
% Date: 10-16-2011
% ========================================================================

function [D,X,T,W]=labelconsistentksvd2(Y,Dinit,Q_train,Tinit,H_train,Winit,iterations,sparsitythres,sqrt_alpha,sqrt_beta)

params.data = [Y;sqrt_alpha*Q_train;sqrt_beta*H_train];
params.Tdata = sparsitythres; % spasity term
params.iternum = iterations;
params.memusage = 'high';
D_ext3 = [Dinit;sqrt_alpha*Tinit;sqrt_beta*Winit];
D_ext3=normcols(D_ext3); % normalization
params.initdict = D_ext3;
% ksvd process
[Dksvd,X,err] = ksvd(params,'');

% get back the desired D, T, W
i_start_D = 1;
i_end_D = size(Dinit,1);
i_start_T = i_end_D+1;
i_end_T = i_end_D+size(Tinit,1);
i_start_W = i_end_T+1;
i_end_W = i_end_T+size(Winit,1);
D = Dksvd(i_start_D:i_end_D,:);
T = Dksvd(i_start_T:i_end_T,:);
W = Dksvd(i_start_W:i_end_W,:);

% normalization
l2norms = sqrt(sum(D.*D,1)+eps);
D = D./repmat(l2norms,size(D,1),1);
T = T./repmat(l2norms,size(T,1),1);
W = W./repmat(l2norms,size(W,1),1);
T = T./sqrt_alpha;
W = W./sqrt_beta;
