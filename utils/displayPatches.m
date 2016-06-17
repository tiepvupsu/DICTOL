function displayPatches2(D, nRows, nCols)
% function tmp = displayPatches(D)
% Display all columns of `D` in an panel 
% require: `size(D, 1)` is a square number or `size(D, 1)/3` is a square number 
% (the former case corresponds to color)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 6/1/2016 10:09:38 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

   V = 1;
   [n, K] = size(D);
   size(D);
   sizeEdge=sqrt(n/V);
   if floor(sizeEdge) ~= sizeEdge
      V=3;
      sizeEdge=sqrt(n/V);
   end
   %% nomalization 
   p=3.5;
   m=min((D(:)));
   if (m >= 0)
      me=0;
      sig=sqrt(mean(((D(:))).^2));
   else
      me=mean(D(:));
      sig=sqrt(mean(((D(:)-me)).^2));
   end
   D = D-me;
   D = min(max(D,-p*sig),p*sig);
   M = max((D(:)));
   m = min((D(:)));
   D = (D-m)/(M-m);
   %% Prepare the panel
%    nRows=floor(sqrt(K)); % number of rows 
%    nCols = ceil(K/nRows); % number of columns
    if nargin == 1 
        [nRows, nCols] = my_factor(K);
    elseif nargin == 2 
        nCols = ceil(K/nRows);
    end
   tmp = zeros((sizeEdge+1)*nRows+1,(sizeEdge+1)*nCols+1,V);
   for r = 1:nRows
      for c = 1:nCols
         if (r-1)*nCols+c > K 
              break;
         end 
         patchCol = D(:,(r-1)*nCols+c);
         patchCol = reshape(patchCol, [sizeEdge,sizeEdge V]);
         tmp(((r-1)*(sizeEdge+1)+2): (r*(sizeEdge+1)),...
             ((c-1)*(sizeEdge+1)+2) : (c*(sizeEdge+1)), :) = patchCol;
      end
   end
   %% display
   colormap('jet');
   imagesc(tmp);
   axis equal off;   
end 

