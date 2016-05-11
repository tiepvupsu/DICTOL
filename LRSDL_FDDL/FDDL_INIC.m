function [opts] = FDDL_INIC (ipts,par)
% ========================================================================
% Coefficient Initialization of FDDL, Version 1.0
% Copyright(c) 2011  Meng YANG, Lei Zhang, Xiangchu Feng and David Zhang
% All Rights Reserved.
%
% ----------------------------------------------------------------------- 
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for initializing the
% Coefficient matrix of FDDL
%
% Please refer to the following paper
%
% Meng Yang, Lei Zhang, Xiangchu Feng, and David Zhang,"Fisher Discrimination 
% Dictionary Learning for Sparse Representation", In IEEE Int. Conf. on
% Computer Vision, 2011.
% L. Rosasco, A. Verri, M. Santoro, S. Mosci, and S. Villa. Iterative
% Projection Methods for Structured Sparsity Regularization. MIT Technical
% Reports, MIT-CSAIL-TR-2009-050,CBCL-282, 2009.
% J. Bioucas-Dias, M. Figueiredo, ?A new TwIST: two-step iterative shrinkage
% /thresholding  algorithms for image restoration?, IEEE Transactions on 
% Image Processing, December 2007.
%----------------------------------------------------------------------
%
%  Inputs :   (1) ipts :    the structre of input data
%                    .D     the dictionary
%                    .X     the training data
%                    .last_coef   the coef in the last iteration
%             (2) par :     the struture of input parameters
%                    .tau   the parameter of sparse constraint of coef
%                    .lambda  the parameter of within-class scatter
%
% Outputs:    (1) opts :    the structure of output data
%                    .A     the coefficient matrix
%                    .ert   the total energy sequence
%
%---------------------------------------------------------------------

par.initM       =     'zero';  % initialization method
par.nIter       =     200;     % maximal iteration number
par.isshow      =     true;    %
par.twist       =     true;    % 'true': use twist
par.citeT       =     1e-6;    %  stop criterion
par.cT          =     1e+10;   %  stop criterion

m    =    size(ipts.D,2);
n    =    size(ipts.X,2);

switch lower(par.initM)
    case {'zero'}
        A    =    zeros(m,n);
    case {'transpose'}
        A  =  ipts.D'*ipts.X;
    case {'pinv'}
        A  =  pinv(ipts.D)*ipts.X;
    case {'last'}
        A    =    ipts.last_coef;
    otherwise
        error('Nonknown method!');
end

D        =    ipts.D;
X        =    ipts.X;
tau      =    par.tau;
lambda   =    par.lambda;
nIter    =    par.nIter;
c        =    par.c;
sigma    =    c;
tau1     =    tau/2;
B        =    eye(n)-ones(n,n)/n;

At_pref   =    A(:);
At_now    =    A(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TWIST parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for_ever           =         1;
IST_iters          =         0;
TwIST_iters        =         0;
sparse             =         1;
verbose            =         1;
enforceMonotone    =         1;
lam1               =         1e-4;   %default minimal eigenvalues
lamN               =         1;      %default maximal eigenvalues
rho0               =         (1-lam1/lamN)/(1+lam1/lamN); 
alpha              =         2/(1+sqrt(1-rho0^2));        %default,user can set
beta               =         alpha*2/(lam1+lamN);         %default,user can set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xm2       =      At_pref;
xm1       =      At_pref;
   
A     =   reshape(At_pref,[m,n]);
gap1  =   norm((X-D*A),'fro')^2;
if n==1
     gap2 = norm(A*B,2)^2;
else
     gap2  =   norm(A*B,'fro')^2;
end
gap3  =   sum(abs(A(:)));
prev_f   =   gap1+2*tau1*gap3+lambda*gap2;
   
for n_it = 2 : nIter;
     A     =   reshape(At_now,[m,n]);
     gap1  =   norm((X-D*A),'fro')^2;
     if n==1
         gap2 = norm(A*B,2)^2;
     else
         gap2  =   norm(A*B,'fro')^2;
     end
     gap3  =   sum(abs(A(:)));
     ert(n_it-1)   =   gap1+2*tau1*gap3+lambda*gap2;
  
%    fprintf('Iteration:%f  Total gap:%f\n',n_it,ert(n_it-1));
    
    while for_ever
     % IPM estimate
        v1   =    [];
        for i  =   1:n
            A     =   reshape(xm1,[m,n]);
            tem1  =   X(:,i)-D*A(:,i);
            tem2  =   D'*tem1;
            v1    =   [v1;tem2];
        end
        A     =   reshape(xm1,[m,n])';
        v2_2   =   [];
        for i  =   1:m
            tem1  =  B*A(:,i);
            v2_2  =  [v2_2;tem1];
        end
        v2_3  =  reshape(v2_2,[n m])';
        v2    =  v2_3(:);
    
        v     =  xm1+(v1-lambda*v2)/sigma;
        x_temp  =  soft(v,tau1/sigma);
        
        if (IST_iters >= 2) | ( TwIST_iters ~= 0)
            % set to zero the past when the present is zero
            % suitable for sparse inducing priors
            if sparse
                mask    =   (x_temp ~= 0);
                xm1     =   xm1.* mask;
                xm2     =   xm2.* mask;
            end
            % two-step iteration
            xm2    =   (alpha-beta)*xm1 + (1-alpha)*xm2 + beta*x_temp;
            % compute residual
            
            A     =   reshape(xm2,[m,n]);
            gap1  =   norm((X-D*A),'fro')^2;
            if n==1
             gap2 = norm(A*B,2)^2;
            else
            gap2  =   norm(A*B,'fro')^2;
            end
            gap3  =   sum(abs(A(:)));
            f   =   gap1+2*tau1*gap3+lambda*gap2;
          
            if (f > prev_f) & (enforceMonotone)
                TwIST_iters   =  0;  % do a IST iteration if monotonocity fails
            else
                TwIST_iters =   TwIST_iters+1; % TwIST iterations
                IST_iters   =    0;
                x_temp      =   xm2;
                if mod(TwIST_iters,10000) ==0
                   c = 0.9*c; 
                   sigma= c;
                end
                break;  % break loop while
            end
        else
          A     =   reshape(x_temp,[m,n]);
          gap1  =   norm((X-D*A),'fro')^2;
          if n==1
            gap2 = norm(A*B,2)^2;
          else
            gap2  =   norm(A*B,'fro')^2;
          end
          gap3  =   sum(abs(A(:)));
          f   =   gap1+2*tau1*gap3+lambda*gap2;

          if f > prev_f
                % if monotonicity  fails here  is  because
                % max eig (A'A) > 1. Thus, we increase our guess
                % of max_svs
                c         =    2*c; 
                sigma     =    c;
                if verbose
%                     fprintf('Incrementing c=%2.2e\n',c);
                end
                if  c > par.cT
                    break;  % break loop while    
                end
                IST_iters = 0;
                TwIST_iters = 0;
           else
                TwIST_iters = TwIST_iters + 1;
                break;  % break loop while
           end
        end
    end

    citerion      =   abs(f-prev_f)/prev_f;
    if citerion < par.citeT | c > par.cT
%        fprintf('Stop!\n c=%2.2e\n citerion=%2.2e\n',c,citerion);
       break;
    end
    
    xm2           =   xm1;
    xm1           =   x_temp;
    At_pref       =   At_now;
    At_now        =   x_temp;
    prev_f        =   f;
    
end

opts.A     =       reshape(At_now,[m,n]);
opts.ert   =       ert;
if par.isshow
    plot(ert,'r-');
end