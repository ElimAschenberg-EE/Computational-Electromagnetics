% ECE 541 Fall 2025 Project 6
% problem 1(c) - Scattered H

clc 
clear

j = 1i;
lambda = 1;
a = 1*lambda;
k = 2*pi/lambda;

M = 40;

x = -5*lambda:0.05*lambda:5*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);

Hz_sc = zeros(2*M+1, length(y), length(x));

ka = k*a;

for n = -M:M

    idx = n + M + 1;

    dJdx = 0.5*( besselj(n-1,ka) - besselj(n+1,ka) );
    dHdx = 0.5*( besselh(n-1,2,ka) - besselh(n+1,2,ka) );

    a_n = - (j^(-n)) * (dJdx / dHdx);

    Hz_sc(idx,:,:) = a_n .* besselh(n,2,k*RHO) .* exp(j*n*PHI);
end

S = sum(Hz_sc, 1);
Hz_sc = squeeze(S(1,:,:));

mask = (RHO < a);

Hz_plot = real(Hz_sc);
Hz_plot(mask) = NaN;

figure;

h = imagesc(x, y, Hz_plot);
caxis([-2 2]);
set(gcf,'Renderer','opengl');                 
set(h, 'AlphaData', ~isnan(Hz_plot));         
set(gca,'Color','w');                   

axis equal tight;
title(sprintf('H_z^{sc}, M = %d', M));
colormap(jet);
colorbar;
