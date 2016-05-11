function [P, Q] = invMyForm(A, B, C)
    % return inv(myForm(A, B, C))
    % -----------------------------------------------
    % Author: Tiep Vu, thv102@psu.edu, 4/8/2016
    %         (http://www.personal.psu.edu/thv102/)
    % -----------------------------------------------
    if nargin == 0
        d = 10;
        C = 10;
        A = rand(d, d);
        B = rand(d, d);
        M = myForm(A, B, C);
    end 
    %%
    invA = inv(A);
    invB = inv(B);
    % Q = myForm(invA, -invA*inv(invB + C*invA)*invA, C);
    P = invA;
    Q = -invA*inv(invB + C*invA)*invA;
    % toc;
    if nargin == 0
        M1 = inv(M); 
        M2 = myForm(P, Q, C);
        norm(M1 - M2)
        pause 
    end 
end 