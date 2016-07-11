function [D, D0, X, X0] = LRSDL_init(Y, Y_range, D_range, opts)
    % opts 
    nClasses = numel(Y_range) - 1;
    D        = zeros(size(Y, 1), D_range(end));
    D0       = zeros(size(Y, 1), opts.k0);    
    X        = zeros(D_range(end), size(Y,2));
    X0       = zeros(opts.k0, size(Y, 2));   
    fprintf('\n');
    for c = 1: nClasses
        if opts.verbose
            fprintf('class %d... ', c);
        end 
        % c
        Yc = get_block_col(Y, c, Y_range);
        switch opts.initmode
            case 'normal'
                [Dc, Xcc] = DL_smallIntraCoeff(Yc, D_range(c+1) - D_range(c), ...
                    opts.lambda1, opts.lambda2, opts);
            otherwise
                Dc = PickDfromY(Yc, [0,size(Yc,2)], ...
                      D_range(c+1) - D_range(c));
                Xcc = pinv(Dc)*Yc;
        end
        %% ========= Xcc = Xcc ==============================        
        col_range = D_range(c) + 1 : D_range(c+1);
        row_range = Y_range(c) + 1 : Y_range(c+1);
        X(col_range, row_range) = Xcc;
        %% ========= D ==============================        
        D(:, col_range) = Dc;
        if opts.verbose
            fprintf('\n');
            if mod(c, 10) == 0
                fprintf('\n');
            end
        end 
    end 
            
    %% ========= Init D0, X0 ==============================    
    if opts.k0 ~= 0
        if opts.verbose
            fprintf('shared dictionary... ');
        end
        Ybar = Y;
        switch opts.initmode
            case 'normal'
                [D0, X0] = LRSDL_initD0X0(Ybar, opts, opts);    
            otherwise
                D0 = PickDfromY(Y, [0, size(Y,2)], opts.k0);
                X0 = pinv(D0)*Y;
        end
    end     
end

