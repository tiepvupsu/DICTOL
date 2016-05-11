    function M = buildMhat(M, M_range_row, M_range_col)
        %% suppose M = [M11 M12 ... M1n;
        %               M21 M22 ... M3n;
        %               .....
        %               Mn1 Mn2 .... Mnn]
        % then Mhat = = [2*M11  M12     ... M1n;
        %                M21    2*M22   ... M3n;
        %                   .....
        %                Mn1     Mn2 .... 2*Mnn]
        C = numel(M_range_row) - 1;
        for c = 1: C
            range_row = M_range_row(c) + 1: M_range_row(c+1);
            range_col = M_range_col(c) + 1: M_range_col(c+1);
            M(range_row, range_col) = 2*M(range_row, range_col);
        end 
    end 