function W = normalizeW(W) % make sum each column of W = 1 
    sw = sum(W, 2);
    W = W./repmat(sw, 1, size(W,2));
end 
