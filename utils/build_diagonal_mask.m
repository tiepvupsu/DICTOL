function a_mask_id = build_diagonal_mask(row_range, col_range)

    C = numel(row_range) - 1;
    M = zeros(row_range(end), col_range(end));
    for c = 1: C 
        M(row_range(c) +1: row_range(c+1), col_range(c)+1: col_range(c+1)) = ...
            ones(row_range(c+1) - row_range(c), col_range(c+1) - col_range(c));
    end 
    a_mask_id = find(M);
end 
