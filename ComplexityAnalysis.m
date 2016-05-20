C = 100;
c = C;
d = 300;
k = 7;
n = 7;
q = 10;

%% old 
c1 = q*k*d^3
%% new 
c2 = n*d*k + n*k^2 + q*(d*k^2 + d^2*k) + d^3 - 3.34e7

c3 = c^2*k*(d*(k+n) + q*n*(c*k))
c4 = c^2*k*n*(d + q*c*k) + C^3*d*k^2
% 
% d = 3000;
% n = 5;
% 
% A = rand(d, d);
% B = rand(d, d);
% 
% tic; 
% for i = 1: n
%     A*B; 
% end  
% toc 
% tic; 
% for i = 1: n
%     inv(A); 
% end
% toc

%% FDDL update D 
% original 
c5 = C*(k*d*n*(C^2 + 1) + k^2*(C*n + q*d))
% efficient 
c6 = C^2*k*(n*d + k*C*n + d*k)

A = C^2*k*n*(d + q*C*k) + C^3*d*k^2
B = C^2*k*d*(k+n) + q*C^2*k*n*(c*k + 1)