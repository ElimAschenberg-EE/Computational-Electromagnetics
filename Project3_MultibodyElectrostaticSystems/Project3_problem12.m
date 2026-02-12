clc
clear

a = 1e-2;
h = 4;
d = 1;
D = sqrt(4*h^2+d^2);
epsilon = 8.854e-12;
sigma = 58e6;
f = 1e6;
mu = 4 * pi * 10^-7;
a11 = 1/(2*pi*epsilon) * log(2*h/a);
a22 = a11;
a12 = 1/(2*pi*epsilon) * log(D/d);
a21 = a12;
aa = a11*a22;
ab = a12*a21;
ac = a11*a22-a12*a21;
C_11 = a22 / (a11*a22-a12*a21);
c = 3e8;

E_cro = 3e6;
Q_cr = E_cro * 2 * pi * epsilon * a;
V_cr = Q_cr/C_11;

w_e = 0.5 * C_11 * V_cr^2;

Rs = sqrt(pi*mu*f/sigma);

R_PUL = Rs / (2 * pi * a)*(1+(-a21/a22)^2);
Z0 = 1/ (c * C_11);
G = sigma / epsilon * C_11;
Y0 = Z0^-1;
alpha = R_PUL/(2*Z0)+G/(2*Y0);
L = epsilon * mu / C_11