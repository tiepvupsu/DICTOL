clear all
fprintf('\nTest myForm function\n');
%%
d = 2;
k = 2;

M = rand(d, d);
N = rand(d, d);
P = rand(d, d);
Q = rand(d, d);
tol = 1e-6;
%%
fprintf('1. Multiplication..............')
A1 = myForm(M, N, k)*myForm(P, Q, k);
[R, S] = myForm_mult(M, N, P, Q, k);
A2 = myForm(R, S, k);
dif = norm(A1 - A2);
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
%% 
fprintf('2. Multiply with a vector......');
x = rand(d*k, 1);
y1 = myForm(M, N, k)*x;
y2 = myForm_mult_vec(M, N, x, k);
dif = norm(y1 - y2);
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
%% 
fprintf('3. Inversion...................');
A1 = inv(myForm(M, N, k));
[R, S] = myForm_inv(M, N, k);
A2 = myForm(R, S, k);
dif = norm(A1 - A2);
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
%% 
fprintf('4. Test LDL problem.......');
F = M*M'; G = N*N';
Q = myForm(F, G, k);
% form big matrix
rho = 10;
e = ones(d,1);
A = myForm(e', zeros(size(e')), k);
QQ = [Q A'/rho; A zeros(size(A,1), size(A,1))];
invQQ = inv(QQ);
% R 
R = inv(Q); % def 
[Fbar, Gbar] = myForm_inv(F, G, k);
R1 = myForm(Fbar, Gbar, k);
fprintf('\n..Check R...');
dif = norm(R - R1);
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
%% S
S = inv(A*R*A');
fbar = sum(Fbar, 2); 
gbar = sum(Gbar, 2);
[S_1, S_2] = myForm_inv(sum(fbar), sum(gbar), k);

fprintf('..Check S...');
dif = norm(S - myForm(S_1, S_2, k));
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
% RAtS
RAtS = R*A'*S;
fprintf('..Check RAtS....');
[RAtS_1, RAtS_2] = myForm_mult(fbar, gbar, S_1, S_2, k);
RAtS2 = myForm(RAtS_1, RAtS_2, k);
dif = norm(RAtS - RAtS2);
if dif < 1e-4
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end
%% RATSAR
RATSAR = RAtS*A*R;
[tmp1, tmp2] = myForm_mult(RAtS_1, RAtS_2, fbar', gbar', k);
RATSAR2 = myForm(tmp1, tmp2, k);
dif = norm(RATSAR - RATSAR2);
fprintf('..Check RAtSAR..');
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 
% N
N = R - RATSAR;
fprintf('..Check N..');
N_1 = Fbar - tmp1;
N_2 = Gbar - tmp2;
N1 = myForm(N_1, N_2, k);
dif = norm(N - N1);
if dif < tol
	fprintf('PASS\n');
else 
	fprintf('FAIL, diff = %5f\n', dif);
end 

invQQ2 = [N RAtS; rho*S*A*R -rho*S];
norm(invQQ - invQQ2)
fprintf('\n');
v = rand(size(Q,1),1);
b = rand(size(A,1),1);
x = [v; b];
z1 = invQQ*x;
z1 = z1(1: size(Q,1));
z2 = myForm_mult_vec(N_1, N_2, v, k) + myForm(RAtS_1, RAtS_2, k)*b;
norm(z1 - z2)
