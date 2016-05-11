% ========================================================================
% Label consistent KSVD algorithm 1
% USAGE: [D,X,T,W]=labelconsistentksvd1(Y,Dinit,Q_train,Tinit,H_train,....
%                        iterations,sparsitythres,sqrt_alpha)
% Inputs
%       Y               -training features
%       Dinit           -initialized dictionary
%       Q_train         -optimal code matrix for training feature 
%       Tinit           -initialized transform matrix
%       H_train         -labels matrix for training feature 
%       iterations      -iterations for KSVD
%       sparsitythres   -sparsity threshold for KSVD
%       sqrt_alpha      -contribution factor
% Outputs
%       D               -learned dictionary
%       X               -sparsed codes
%       T               -learned transform matrix
%       W               -learned classifier parameters
%
% Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
% Date: 10-16-2011
% ========================================================================

function [D,X,T,W]=labelconsistentksvd1(Y,Dinit,Q_train,Tinit,H_train,iterations,sparsitythres,sqrt_alpha)

params.data = [Y;sqrt_alpha*Q_train];
params.Tdata = sparsitythres; % spasity term
params.iternum = iterations;
params.memusage = 'high';
D_ext2 = [Dinit;sqrt_alpha*Tinit];
D_ext2=normcols(D_ext2); % normalization
params.initdict = D_ext2;
% ksvd process
[Dksvd,X,err] = ksvd(params,'');

% get back the desired D, T
i_start_D = 1;
i_end_D = size(Dinit,1);
i_start_T = i_end_D+1;
i_end_T = i_end_D+size(Tinit,1);
D = Dksvd(i_start_D:i_end_D,:);
T = Dksvd(i_start_T:i_end_T,:);

% normalization
l2norms = sqrt(sum(D.*D,1)+eps);
D = D./repmat(l2norms,size(D,1),1);
T = T./repmat(l2norms,size(T,1),1);
T = T./sqrt_alpha;

% learning linear classifier parameters
W = inv(X*X'+eye(size(X*X')))*X*H_train';
W = W';