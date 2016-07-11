function new_range = range_delete_ids(a_range, ids)
% function new_range = range_delete_ids(a_range, ids)
% given a range `a_range` of an array. Suppose we want to delete some
% element of that array indexed by `ids`, `new_range` is the new range


    if nargin == 0 %% test input
        a_range = [0, 3, 5, 10];
        ids = [1, 4, 2, 7, 10];
    end 
    %% MAIN
    ids = sort(ids);
    n = numel(a_range);
    a = zeros(1, n);
    j = 2;
    while j < n
        for i = 1: numel(ids)
            while a_range(j) < ids(i)
                j = j + 1;
            end 
            for k = j:n
                a(k) = a(k) + 1;
            end 
        end 
    end 
    new_range = a_range - a;
    %% test output
    if nargin == 0
        disp(a_range);
        disp(a);
        disp(new_range);
        new_range = [];
    end 
end 