% ========================================================================
% Classification 
% USAGE: [prediction, accuracy, err] = classification(D, W, data, Hlabel,
%                                       sparsity)
% Inputs
%       D               -learned dictionary
%       W               -learned classifier parameters
%       data            -testing features
%       Hlabel          -labels matrix for testing feature 
%       iterations      -iterations for KSVD
%       sparsity        -sparsity threshold
% outputs
%       prediction      -predicted labels for testing features
%       accuracy        -classification accuracy
%       err             -misclassfication information 
%                       [errid featureid groundtruth-label predicted-label]
%
% Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
% Date: 10-16-2011
% ========================================================================

function [prediction, accuracy, err] = classification(D, W, data, Hlabel, sparsity)

% sparse coding
G = D'*D;
% Gamma = omp(D'*data,G,sparsity);
Gamma = myOMP(data, D, sparsity);
% classify process
errnum = 0;
err = [];
prediction = [];
for featureid=1:size(data,2)
    spcode = Gamma(:,featureid);
    score_est =  W * spcode;
    score_gt = Hlabel(:,featureid);
    [maxv_est, maxind_est] = max(score_est);  % classifying
    [maxv_gt, maxind_gt] = max(score_gt);
    prediction = [prediction maxind_est];
    if(maxind_est~=maxind_gt)
        errnum = errnum + 1;
        err = [err;errnum featureid maxind_gt maxind_est];
    end
end
accuracy = (size(data,2)-errnum)/size(data,2);