function H = lcksvd_buildH(label)
    Y_range = label_to_range(label);
    C = numel(Y_range) - 1;
    H = zeros(C, numel(label));

    for c = 1: C 
        H(c, Y_range(c) + 1: Y_range(c+1)) = ones(1, Y_range(c+1) - Y_range(c));
    end 
end 
