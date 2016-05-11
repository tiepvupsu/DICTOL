function D    =    FDDL_INID(data,nCol,wayInit)
% ========================================================================
% Dictionary Initialization of FDDL, Version 1.0
% Copyright(c) 2011  Meng YANG, Lei Zhang, Xiangchu Feng and David Zhang
% All Rights Reserved.
%
% -----------------------------------------------------------------------    
%   
% Input :   (1) data :  the data matrix 
%           (2) nCol :  the number of dictioanry's columns
%           (3) wayInit:  the method to initialize dictionary
% 
% Outputs : (1) D  :    the initialized dictionary
%
%------------------------------------------------------------------------

m   =    size(data,1);

switch lower(wayInit)
    case {'pca'}
        [D,disc_value,Mean_Image]   =    Eigenface_f(data,nCol-1);
        D                           =    [D Mean_Image./norm(Mean_Image)];
    case {'random'}
        phi                         =    randn(m, nCol);
        phinorm                     =    sqrt(sum(phi.*phi, 2));
        D                           =    phi ./ repmat(phinorm, 1, nCol);
    otherwise 
        error{'Unkonw method.'}
end
return;
