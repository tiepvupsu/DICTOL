function cost = DLSI_term(D, D_range)
% * Syntax: cost = DLSI_term(D, D_range)
% * Calculating the structured incoherence term in DLSI [[5]](#fn_dls).
% * $\sum_{c=1}^C \sum_{i \neq c} \|D_i^TD_c\|_F^2$
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
        d = 100;
        D_range = 10*(0:10);
        D = rand(d, D_range(end));
    end 
    %% MAIN
    A = erase_diagonal_blocks(D'*D, D_range, D_range);
    cost = normF2(A);
    %% 
    if nargin == 0 
        cost = [];
    end 
end 