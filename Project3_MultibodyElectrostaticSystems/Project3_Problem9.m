clc
clear

d = 2e-2;
S = 1;
epsilon = 8.854e-12;

a11 = 3*d/(epsilon * S);
a12 = 2*d/(epsilon * S);
a13 = d/(epsilon * S);
a21 = a12;
a22 = 2*d/(epsilon * S);
a23 = d/(epsilon * S);
a31 = a13;
a32 = a23;
a33 = d/(epsilon * S);

a = [a11, a12, a13; 
     a21, a22, a23; 
     a31, a32, a33];
b_matrix = inv(a);

c11 = b_matrix(1,1)+b_matrix(1,2) + b_matrix(1,3);
c22 = b_matrix(2,1) + b_matrix(2,2) + b_matrix(2,3);
c33 = b_matrix(3,1) + b_matrix(3,2) + b_matrix(3,3);
c12 = -b_matrix(1,2);
c13 = -b_matrix(3,1);
c21 = -b_matrix(2,1);
c23 = -b_matrix(2,3);
c31 = -b_matrix(3,1);
c32 = -b_matrix(3,2);
c = [c11 c12 c13;
     c21 c22 c23;
     c31 c32 c33];
Q = 10e-6;
V = 10e3;
E_12 = Q/(epsilon * S);
E_13 = Q/(epsilon * S);
E_23 = Q/(epsilon * S);