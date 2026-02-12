% ECE 541 Fall 2025 Project 6
% problem 1(c) - Total H

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

Hz_tot = zeros(size(RHO));

for n = -M:M

    J_ka_prime = 0.5*( besselj(n-1,ka) - besselj(n+1,ka) );
    H_ka_prime = 0.5*( besselh(n-1,2,ka) - besselh(n+1,2,ka) );

    b_n = -j^(-n) * (J_ka_prime / H_ka_prime);

    Hz_inc_n = j^(-n) .* besselj(n, k*RHO) .* exp(j*n*PHI);
    Hz_sc_n  = b_n .* besselh(n,2,k*RHO) .* exp(j*n*PHI);

    Hz_tot = Hz_tot + (Hz_inc_n + Hz_sc_n);
end

mask = (RHO < a);
Hz_plot = real(Hz_tot);
Hz_plot(mask) = NaN;

figure;
h = imagesc(x, y, Hz_plot);
set(gcf,'Renderer','opengl');            
set(h,'AlphaData', ~isnan(Hz_plot));    
set(gca,'Color','w');                

axis equal tight;
title(sprintf('Re{H_z^{tot}} (TE_z), M = %d', M));
colormap(jet);
colorbar;
caxis([-2 2]);
