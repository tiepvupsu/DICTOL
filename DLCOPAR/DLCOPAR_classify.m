function pred = DLCOPAR_classify(Y, D, opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: test mode ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
        d = 100;
        C = 10;
        N = 10;
        k = 10;
        k0 = 20;
        Y = normc(rand(d, C*N));
        D = normc(rand(d, C*k + k0));
        opts.D_range_ext = [k*(0:C) C*k+k0];
        if rand(1) < 0.5
            opts.classify_mode = 'LC';
        else 
            opts.classify_mode = 'GC';
        end 
        opts.gamma = 0.01;
    end         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: test mode ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    if isfield(opts, 'classify_mode') == 0
        fprintf('You need to specify classification mode in opts.classify_mode: GC/LC\n');
    elseif strcmp(opts.classify_mode, 'GC')
        pred = GC(Y, D, opts);
    elseif strcmp(opts.classify_mode, 'LC')
        pred = LC(Y, D, opts);
    else 
        fprintf('classify_mode is either GC or LC');
    end             
end 

function pred = GC(Y, D, opts)
%     fprintf('GC mode\n');
    C = numel(opts.D_range_ext) - 2;
    D_range_ext = opts.D_range_ext;
    pars.show = false;
    pars.max_iter = 100;
    if isfield(opts, 'gamma') == 0
        fprintf('specify gamma');
    else 
%         X = myLasso_fista_2(Y, D, zeros(size(D,2), size(Y, 2)), opts.gamma, pars);
        X = lasso_fista(Y, D, zeros(size(D, 2), size(Y, 2)), opts.gamma, pars);
        DCp1_range = D_range_ext(C+1)+1: D_range_ext(C+2); 
        DCp1 = D(:, DCp1_range);
        XCp1 = X(DCp1_range, :);
        Y = Y - DCp1*XCp1;
        E = zeros(C, size(Y,2));
        for i = 1:C
            Xi = get_block_row(X, i, D_range_ext);
            Di = get_block_col(D, i, D_range_ext);
            R = Y - Di*Xi;
            E(i,:) = sum(R.^2, 1);
        end
        [~,pred] = min(E);
    end 
    
end 


function pred = LC(Y, D, opts)
%     fprintf('LC mode\n');
    C = numel(opts.D_range_ext) - 2;
    D_range_ext = opts.D_range_ext;
    pars.show = false;
    pars.max_iter = 100;
    if isfield(opts, 'gamma') == 0
        fprintf('specify gamma');
    else 
        DCp1_range = D_range_ext(C+1)+1: D_range_ext(C+2); 
        E = zeros(C, size(Y,2));
        for i = 1:C
            Dc_range = D_range_ext(i)+1: D_range_ext(i+1);
            Dchat = D(:, union(Dc_range, DCp1_range));
%             X = myLasso_fista_2(Y, Dchat, zeros(size(Dchat,2), size(Y, 2)), opts.gamma, pars);
            X = lasso_fista(Y, Dchat, zeros(size(Dchat,2), size(Y, 2)), opts.gamma, pars);

            R1 = Y - Dchat*X;
            R2 = opts.gamma*abs(X);
            E(i,:) = sum(R1.^2, 1) + sum(R2.^2, 1);
        end
        [~,pred] = min(E);
    end 
end 