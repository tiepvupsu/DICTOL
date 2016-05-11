function X = DLCOPAR_updateX(Y, Y_range, D, X, opts) %page 189-190 DLCOPAR
    % see DLCOPAR paper: http://www.cs.zju.edu.cn/people/wangdh/papers/draft_ECCV12_particularity.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ================== block: test module 12/19/2015 4:03:53 PM =========================%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    fprintf('Test most\n');
%     profile off;
%     profile on;
    d = 3000;
    N = 30;
    k = 20;
    k0 = 50;
    C = 100 ;
    opts.k0 = k0;
    opts.lambda = 0.01;
    opts.eta = 0.1;
    Y = normc(rand(d, C*N));
    D = normc(rand(d, C*k + k0));
    X = zeros(size(D,2), size(Y,2));
%     train_label = [];
%     for c = 1: C
%         train_label = [train_label c*ones(1, N)];
%     end         
    Y_range = N*(0:C);
    opts.D_range = k* (0:C);
    opts.D_range_ext = [opts.D_range opts.D_range(end)+k0];
    opts.show = false;
    opts.max_iter = 30;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ------------------end of block: test module 12/19/2015 4:03:53 PM -------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    C = numel(Y_range) - 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ================== block: fista aproach ==========================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DtD = D'*D;
    DtY = D'*Y;
    L = 2*max(eig(DtD)) + 10;
    
    optsX = opts;
    optsX.max_iter = 100;

    for c = 1: C
        % fprintf('class: %3d | ', c);
        Xc = get_block_col(X, c, Y_range);
        X(:, Y_range(c)+1: Y_range(c+1)) = DLCOPAR_updateXc(DtD, DtY, Y_range, Xc, c, L, optsX);
        if nargin == 0 % test mode
            fprintf('class = %3d || cost: %5f\n', c, DLCOPAR_cost(Y, Y_range, D, X, opts));
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% ------------------end of block: fista aproach ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%     fprintf('after cost: %5f\n', DLCOPAR_cost(Y, Y_range, D, X, opts));

    %% ========= view profile - for testing ==============================    
    if nargin == 0
        profile viewer;
        return;
    end
end 

%%
% function Qcom_c_tilder = buildQcom_c_tilder(c, D_range_ext) % page 189 DLCOPAR
%     C = numel(D_range_ext) - 2;
%     Qcom_c_tilder = eye(D_range_ext(end));
%     Qcom_c_tilder(:, D_range_ext(C+1)+1: D_range_ext(C+2)) = [];
%     Qcom_c_tilder(:, D_range_ext(c)+1: D_range_ext(c+1)) = [];
% end 

%5
% function Dprime = buildDprime(D, c, D_range_ext) % page 190, DLCOPAR 
%     d = size(D,1);
%     C = numel(D_range_ext) - 2;
%     Qcom_c_tilder = buildQcom_c_tilder(c, D_range_ext);
%     Dprime = [D; zeros(size(D)); Qcom_c_tilder'];
% 
%     Dprime(d+1:2*d, D_range_ext(c)+1: D_range_ext(c+1)) = D(:, D_range_ext(c)+1: D_range_ext(c+1));
%     Dprime(d+1:2*d, D_range_ext(C+1)+1: D_range_ext(C+2)) = D(:, D_range_ext(C+1)+1: D_range_ext(C+2));
% end 