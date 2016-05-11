function [A, B] = mult_myForm(M, N, P, Q, C)
    % return myForm(M, N, C) * myForm(P, Q, C) = myForm(A, B, C);
    % -----------------------------------------------
    % Author: Tiep Vu, thv102@psu.edu, 4/8/2016
    %         (http://www.personal.psu.edu/thv102/)
    % -----------------------------------------------
    if nargin == 0
        d = 10;
        C = 10;
        M = rand(d, d);
        N = rand(d, d);
        P = rand(d, d);
        Q = rand(d, d);
        
    end 
    %%
    A = M*P;
    B = M*Q + N*P + C*N*Q;
    % From here, we can see that if N == 0 => A = M*P; B = M*Q;
    % If Q == 0 => A = M*P; B = N*P.
    % toc;
    if nargin == 0
        M1 = myForm(A, B, C);
        M2 = myForm(M, N, C)*myForm(P, Q, C);
        norm(M1 - M2)
        pause 
    end 
end 