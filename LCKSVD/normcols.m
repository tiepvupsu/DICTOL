function matout=normcols(matin)
l2norms = sqrt(sum(matin.*matin,1)+eps);
matout = matin./repmat(l2norms,size(matin,1),1);