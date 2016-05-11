function X = myLassoWIntrasmall_fista(Y, D, lambda1, lambda2, Xinit, pars)
%% cost = .5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildM(X)) + lambda1*norm1(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ================== block: Test module ==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function myLassoWIntrasmall_fista()
%     d = 30;
%     k = 7;
%     N = 7;
%     Y = normc(rand(d, N));
%     D = normc(rand(d, k));
%     Y = D;

%     lambda1 = 0.01;
%     lambda2 = 1;

%     Xinit = zeros(size(D,2), size(Y, 2));   
%     pars.max_iter = 100;
%     pars.show = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ------------------end of block: Test module ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

	%% =========  ==============================
    
	DtD = D'*D;
	DtY = D'*Y;
	I = eye(size(D,2));
	D1 = DtD + lambda2*I;
	function M = buildM(X)
		M = repmat(mean(X,2), 1, size(X,2));
	end
	function cost = cal_F(X)
		cost = .5*normF2(Y - D*X) + .5*lambda2*normF2(X - buildM(X)) + lambda1*norm1(X);
	end 

	function grad = cal_grad_f(X)
		grad = D1*X - DtY - lambda2*buildM(X);
	end 

	
	k = 0;
	
    % L = max(eig(DtD)) + 1;
	L = max(eig(D1));
	cost_old = 0;
	tol = 1e-4;
    t = 1/L;
    x_init = Xinit;
	x_old = x_init;
    y_old = x_init;
    t_old = 1;
	while k < pars.max_iter
        k = k + 1;
        x_new = shrinkage(y_old - cal_grad_f(y_old)/L, lambda1/L);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
		cost_new = cal_F(x_new);
        % if pars.show
        %     disp(cost_new);
        % end
%         abs(cost_new - cost_old)
		% if(abs(cost_new - cost_old) < tol)
		% 	break;
		% end
        if norm1(x_new - x_old)/numel(x_new) < tol 
            break;
        end

		% cost_old = cost_new;
		x_old = x_new;
		t_old = t_new;
		y_old = y_new;
    end
    X = x_new;
end
