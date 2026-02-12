clc
clear

a = 1e-2;
H = 7;
h = 3;
epsilon = 8.854e-12;

a11 = 1/(2*pi*epsilon) * log(2*H/a);
a12 = 1/(2*pi*epsilon) * log((H+h)/(H-h));
a21 = a12;
a22 = 1/(2*pi*epsilon) * log(2*h/a);
a_matrix = [a11 a12; a21 a22];
b_matrix = a_matrix^-1;
c11 = b_matrix(1,1) + b_matrix(2,1);
c12 = -b_matrix(1,2);
c21 = -b_matrix(2,1);
c22 = b_matrix(2,1) + b_matrix (2,2);
c_matrix = [c11 c12; c21 c22];
aQ1_dif = a11-a21;
aQ2_dif = a22 - a12;
frac1 = aQ1_dif/aQ2_dif;
frac2 = aQ2_dif/aQ1_dif;
C = 1/(a11*frac1 + a12) + 1/(a11+ frac2 * a12)
