function [D, X] = COPAR_init(Y, Y_range, opts)    
% function [D, X] = COPAR_init(Y, Y_range, opts)    
% COPAR initialization
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    C           = numel(Y_range) - 1;
    D_range     = opts.k*(0:C);
    D_range_ext = [D_range D_range(end)+opts.k0];    
    D           = zeros(size(Y,1), D_range_ext(end));
    X           = zeros(D_range_ext(end), Y_range(end));    
    if opts.verbose 
        fprintf('Initializing...\n');
        fprintf('class: \n');
    end 
    for c = 1: C 
        if opts.verbose
            fprintf('%3d ', c);
            if mod(c, 10) == 0
                fprintf('\n');
            end 
        end 
        Yc        = get_block_col(Y, c, Y_range);
        [Dc, Xcc] = ODL(Yc, D_range(c+1) - D_range(c), opts.lambda, opts);
        D(:, D_range(c)+1: D_range(c+1)) = Dc;
        X(D_range(c)+1: D_range(c+1), Y_range(c)+1: Y_range(c+1)) = Xcc;
    end 
    if opts.k0 > 0
        [DCp1, XCp1]                  = ODL(Y, opts.k0, opts.lambda, opts);
        D(:, D_range_ext(C+1)+1: end) = DCp1;
        X(D_range_ext(C+1)+1:end, :)  = XCp1;
    end 
        
end 