% ========================================================================
% Initialization for Label consistent KSVD algorithm
% USAGE: [Dinit,Tinit,Winit,Q] = initialization4LCKSVD(training_feats,....
%                               H_train,dictsize,iterations,sparsitythres)
% Inputs
%       training_feats  -training features
%       H_train         -label matrix for training feature 
%       dictsize        -number of dictionary items
%       iterations      -iterations
%       sparsitythres   -sparsity threshold
% Outputs
%       Dinit           -initialized dictionary
%       Tinit           -initialized linear transform matrix
%       Winit           -initialized classifier parameters
%       Q               -optimal code matrix for training features 
%   
% Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
% Date: 10-16-2011
% ========================================================================


function [Dinit,Tinit,Winit,Q]=initialization4LCKSVD(training_feats,H_train,dictsize,iterations,sparsitythres)

numClass = size(H_train,1); % number of objects
numPerClass = round(dictsize/numClass); % initial points from each classes
Dinit = []; % for LC-Ksvd1 and LC-Ksvd2
dictLabel = [];
for classid=1:numClass
    col_ids = find(H_train(classid,:)==1);
    data_ids = find(colnorms_squared_new(training_feats(:,col_ids)) > 1e-6);   % ensure no zero data elements are chosen
    %    perm = randperm(length(data_ids));
    perm = [1:length(data_ids)]; 
    %%%  Initilization for LC-KSVD (perform KSVD in each class)
    Dpart = training_feats(:,col_ids(data_ids(perm(1:numPerClass))));
    para.data = training_feats(:,col_ids(data_ids));
    para.Tdata = sparsitythres;
    para.iternum = iterations;
    para.memusage = 'high';
    % normalization
    para.initdict = normcols(Dpart);
    % ksvd process
    [Dpart,Xpart,Errpart] = ksvd(para,'');
    Dinit = [Dinit Dpart];
    
    labelvector = zeros(numClass,1);
    labelvector(classid) = 1;
    dictLabel = [dictLabel repmat(labelvector,1,numPerClass)];
end

% Q (label-constraints code); T: scale factor
T = eye(dictsize,dictsize); % scale factor
Q = zeros(dictsize,size(training_feats,2)); % energy matrix
for frameid=1:size(training_feats,2)
    label_training = H_train(:,frameid);
    [maxv1,maxid1] = max(label_training);
    for itemid=1:size(Dinit,2)
        label_item = dictLabel(:,itemid);
        [maxv2,maxid2] = max(label_item);
        if(maxid1==maxid2)
            Q(itemid,frameid) = 1;
        else
            Q(itemid,frameid) = 0;
        end
    end
end

params.data = training_feats;
params.Tdata = sparsitythres; % spasity term
params.iternum = iterations;
params.memusage = 'high';

% normalization
params.initdict = normcols(Dinit);

% ksvd process
[Dtemp,Xtemp,Errtemp] = ksvd(params,'');

% learning linear classifier parameters
Winit = inv(Xtemp*Xtemp'+eye(size(Xtemp*Xtemp')))*Xtemp*H_train';
Winit = Winit';

Tinit = inv(Xtemp*Xtemp'+eye(size(Xtemp*Xtemp')))*Xtemp*Q';
Tinit = Tinit';
