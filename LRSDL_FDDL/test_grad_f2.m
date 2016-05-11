function test_grad_f2()
    addpath(fullfile('..', 'utils'));
    d = 50;
    Y_range = [0, 5, 8, 9, 10];
    X = rand(d, Y_range(end));
    C = numel(Y_range) - 1;
    %%
    function Mc = buildM(Xc)
        Mc = repmat(mean(Xc, 2), 1, size(Xc, 2));
    end 
    %%
    function M = buildMtotal(X)
        M = zeros(size(X));
        for c = 1: C 
            range_c = get_range(Y_range, c);
            M(:, range_c) = buildM(get_block_col(X, c, Y_range));
        end 
    end 
    %%
    function cost = cacl_f(X)
        cost = normF2(X);
        m = mean(X, 2);
        for c = 1: C 
            Xc = get_block_col(X, c, Y_range);
            Mc = buildM(Xc);
            M = repmat(m, 1, size(Mc, 2));
            cost = cost + normF2(Xc - Mc) - normF2(Mc - M);
        end 
        cost = 0.5*cost;

    end 
    %%
    function g = grad(X)
        g = 2*X + buildM(X) - 2*buildMtotal(X);
    end 
    %% 
    check_grad(@cacl_f, @grad, X);
end 
