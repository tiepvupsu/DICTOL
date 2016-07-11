function [X, iter] = fista(grad, Xinit, L, lambda, opts, calc_F)   
% function [X, iter] = fista(grad, Xinit, L, lambda, opts, calc_F)   
% * A Fast Iterative Shrinkage-Thresholding Algorithm for 
% Linear Inverse Problems.
% * Solve the problem: `X = arg min_X F(X) = f(X) + lambda||X||_1` where:
%   - `X`: variable, can be a matrix.
%   - `f(X)` is a smooth convex function with continuously differentiable 
%       with Lipschitz continuous gradient `L(f)` (Lipschitz constant of 
%       the gradient of `f`).
% * Syntax: `[X, iter] = FISTA(calc_F, grad, Xinit, L, lambda,  opts)` where:
%   - INPUT:
%     + `grad`: a _function_ calculating gradient of `f(X)` given `X`.
%     + `Xinit`: initial guess.
%     + `L`: the Lipschitz constant of the gradient of `f(X)`.
%     + `lambda`: a regularization parameter, can be either a scalar or 
%       a weighted matrix.
%     + `opts`: a _structure_ variable describing the algorithm.
%       * `opts.max_iter`: maximum iterations of the algorithm. 
%           Default `300`.
%       * `opts.tol`: a tolerance, the algorithm will stop if difference 
%           between two successive `X` is smaller than this value. 
%           Default `1e-8`.
%       * `opts.show_progress`: showing `F(X)` after each iteration or not. 
%           Default `false`. 
%     + `calc_F`: optional, a _function_ calculating value of `F` at `X` 
%       via `feval(calc_F, X)`. 
% -------------------------------------
% Author: Tiep Vu, thv102, 4/6/2016
% (http://www.personal.psu.edu/thv102/)
% -------------------------------------
    opts = initOpts(opts);
    Linv = 1/L;    
    lambdaLiv = lambda*Linv;
    x_old = Xinit;
    y_old = Xinit;
    t_old = 1;
    iter = 0;
    cost_old = 1e10;
    %% MAIN LOOP
    while  iter < opts.max_iter
        iter = iter + 1;
        x_new = shrinkage(y_old - Linv*feval(grad, y_old), lambdaLiv);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
        %% check stop criteria
        e = norm1(x_new - x_old)/numel(x_new);
        if e < opts.tol
            break;
        end
        %% update
        x_old = x_new;
        t_old = t_new;
        y_old = y_new;
        %% show progress
        if opts.verbose
            if nargin ~= 0
                cost_new = feval(calc_F, x_new);
%                 if cost_new <= cost_old 
%                     stt = 'YES.';
%                 else 
%                     stt = 'NO, check your code.';
%                 end
%                 fprintf('iter = %3d, cost = %f, cost decreases? %s\n', ...
%                     iter, cost_new, stt);
                fprintf('iter = %3d, cost = %f\n', iter, cost_new);
                cost_old = cost_new;
            else 
                if mod(iter, 5) == 0
                    fprintf('.');
                end
                if mod(iter, 10) == 0 
                   fprintf('%d', iter);
                end     
            end        
        end 
    end
    X = x_new;
end 