clc
clear

b_matrix_ansys=[118.022033 -6.156154;
                -6.156154 117.764573]*10^-12;
L_matrix_ansys=[292.681725 35.968671;
                35.968671 292.681243]*10^-9;

a_matrix_ansys = L_matrix_ansys * b_matrix_ansys;

a11 = a_matrix_ansys(1,1);
a22 = a_matrix_ansys(2,2);
a_sum = a11+a22;

syms c_m

eqn = (1/(c_m^2))^2 - (a11+a22)*(1/(c_m^2))+a11*a22 == 0;
cm_sol = solve(eqn, c_m);
cm_num = double(cm_sol)