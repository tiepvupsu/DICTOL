function M = inv_IpXY(X, Y)
% Calculate the inverse of matrix A = I + XY.
% if X is a fat matrix (number of columns > number of rows), then use inv(I + X*Y)
% else: use equation: (I + XY)^(-1) = I - Y*(I + Y*X)^(-1)*X 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/12/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    if nargin == 0
        d1 = 3000;
        d2 = 30;
        X = rand(d1, d2);
        Y = rand(d2, d1);
        tic 
    end 
    %%
    [d1, d2] = size(X);
    if d1 > d2 
        M = inv(eye(d1) + X*Y);
    else 
        M = eye(d1) - X*inv(eye(d2) + Y*X)*Y;
    end
    %% test if no input
    if nargin == 0
        toc 
        tic 
        M1 = inv(eye(d1) + X*Y);
        toc
        norm(M1 - M) % should be close to zero.
        pause 
    end

end 
