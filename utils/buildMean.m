function M = buildMean(X)
    N = size(X, 2);
    m = mean(X, 2);
    M = repmat(m, 1, N);
end 