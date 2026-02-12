clc
clear

a = 1e-3;
b = 2e-3;
c = 2.5e-3;
d = 5e-3;
epsilon0 = 8.854e-12;
epsilonr = 4;
epsilon = epsilonr * epsilon0;
V1 = 3;
V2 = 1;

a11 = (log(b/a)+log(d/c))/(2*pi*epsilon);
a21 = (log(d/c))/(2*pi*epsilon);
a12 = a21;
a22 = (log(d/c))/(2*pi*epsilon);

a_matrix = [a11 a12 ; a21 a22];
b_matrix = inv(a_matrix);
b11 = b_matrix(1,1);
b12 = b_matrix(1,2);
b21 = b_matrix(2,1);
b22 = b_matrix(2,2);
c11 = b_matrix(1,1)+b_matrix(1,2);
c22 = b_matrix(2,1) + b_matrix(2,2);
c12 = -b_matrix(1,2);
c21 = -b_matrix(2,1);
c_matrix = [c11 c12; c21 c22];

Q1 = b11 * V1 + b12 * V2;
Q2 = b21 * V1 + b22 * V2;

We = (1/2) * (Q1*V1+Q2*V2);