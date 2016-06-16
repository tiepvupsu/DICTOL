d = 500;
n = 20;
C = 100;
k = 10;
q = 50;
q2 = 50;

ODLSI = C*k*(k*d+d*n + q*k*n) + C*q*k*d^3
EDLSI = C*k*(k*d+d*n + q*k*n) + C*d^3 + C*q*d*k*(q*k + d)
OFDDL = C^2*d*k*(n+C*k+C*n)+C*k^2*q*(d + C^2*n)
EFDDL = C^2*k*((q+1)*k*(d+C*n) + 2*d*n)

OCOPAR = C^3*k^2*(2*d + C*k + q*n) + C*q*k*d^3
ECOPAR = C^3*k^2*(2*d + C*k + q*n) + C*d^3 + C*q*d*k*(q*k + d)

LRSDL = C^2*k*((q+1)*k*(d + C*n) + 2*d*n) + C^2*d*k*n + (q + q2)*d*k^2

%% 
ODLSID = C*q*k*d^3 
EDLSID = C*d^3 + C*q*d*k*(q*k + k)
OFDDLX = C^2*k*(d*n + q*C*k*n + C*d*k)
EFDDLX = C^2*k*(d*n + q*C *n*k + d*k)
OFDDLD = C*d*k*(q*k + C^2*n)
EFDDLD = C*d*k*(C*n + C*q*k) + C^3*k^2 *n