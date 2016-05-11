function [D, D0, X, X0] = LRSDL_init2(Y, Y_range, D_range, opts, pars )

    D = [];
    X = zeros(D_range(end), size(Y,2));
    nClasses = numel(Y_range) - 1;
    % parstmp.max_iter = 30;    parstmp.show = false;
    parstmp = pars;
    % for c = 1: nClasses
    %     if pars.show
    %         fprintf('Class %3d: ', c);
    %     end 
    %     c
    %     Yc = get_block_col(Y, c, Y_range);
    %     switch opts.initmode
    %         case 'normal'
    %             [Dc, Xcc] = DL_smallIntraCoeff(Yc, D_range(c+1) - D_range(c), opts.lambda1, opts.lambda2, parstmp);
    %         % case 'pca'
    %         %     [Dc, Xcc] = Dict_PCA_init(Yc, opts.k);
    %         otherwise
    %             Dc = PickDfromY(Yc, [0,size(Yc,2)], opts.k);
    %             Xcc = pinv(Dc)*Yc;
    %     end
    %     % [Dc, Xcc] = SDDL_Dec16_init(Yc, D_range(c+1) - D_range(c), opts.lambda2, pars);

    %     %% ========= Xcc = Xcc ==============================        
    %     col_range = D_range(c) + 1 : D_range(c+1);
    %     row_range = Y_range(c) + 1 : Y_range(c+1);
    %     X(col_range, row_range) = Xcc;

    %     %% ========= D ==============================        
    %     D = [D Dc];
    %     if pars.show 
    %         fprintf('\n');
    %     end 
    % end 

    [D, X] = FDDL_init2(Y, Y_range, opts);
    % --------------- Init D0, X0 -------------------------
    % [D0, X0] = DL_smallIntraCoeff(Y, opts.k0, 0.001, opts.lambda2, max_iter);
    % Ybar = Y - D*X;
    Ybar = Y;
    switch opts.initmode
        case 'normal'
            [D0, X0] = LRSDL_initD0X0(Ybar, opts, pars);    
        % case 'pca'
        %     [D0, X0] = Dict_PCA_init(Y, opts.k0);
        otherwise
            D0 = PickDfromY(Y, [0, size(Y,2)], opts.k0);
            X0 = pinv(D0)*Y;
    end
    
    % [D0, X0] = SDDL_Dec16_init(Y, opts.k0, opts.lambda2, pars);
end

