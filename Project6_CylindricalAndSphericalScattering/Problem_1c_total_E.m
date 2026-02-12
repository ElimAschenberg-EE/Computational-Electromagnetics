% ECE 541 Fall 2025 Project 6
% problem 1(c) - Total E

clc
clear

j = 1i;

lambda = 1;
a = 1*lambda;
k = 2*pi/lambda;
ka = k*a;

M = 40;

x = -5*lambda:0.05*lambda:5*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);

Ez_tot = zeros(size(RHO));

for n = -M:M
    c_n = j^(-n);

    b_n = -c_n * ( besselj(n,ka) / besselh(n,2,ka) );

    Ez_inc_n = c_n .* besselj(n, k*RHO) .* exp(j*n*PHI);
    Ez_sc_n  = b_n .* besselh(n, 2, k*RHO) .* exp(j*n*PHI);

    Ez_tot = Ez_tot + (Ez_inc_n + Ez_sc_n);
end

figure;
imagesc(x, y, real(Ez_tot));
axis equal tight;
caxis([-2 2]);
title(sprintf('E_z^{tot} (TM_z), M = %d', M));
colormap(jet)
colorbar;

hold on
th = linspace(0,2*pi,500);
fill(a*cos(th), a*sin(th), 'w', 'EdgeColor','none');
hold off
