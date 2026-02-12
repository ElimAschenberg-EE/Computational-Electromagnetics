% ECE 541 Fall 2025 Project 6
% problem 1(c) - Scattered E

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

Ez_sc = zeros(2*M+1, length(y), length(x));

ka = k*a;

for n = -M:M
    idx = n + M + 1;

    Jn_ka = besselj(n, ka);
    Hn_ka = besselh(n, 2, ka);

    a_n = - (j^(-n)) * (Jn_ka / Hn_ka);

    Ez_sc(idx,:,:) = a_n .* besselh(n, 2, k*RHO) .* exp(j*n*PHI);
end

S = sum(Ez_sc, 1);
Ez_sc = squeeze(S(1,:,:));


mask = (RHO < a);
Ez_plot = real(Ez_sc);
Ez_plot(mask) = NaN;

figure;
h = imagesc(x, y, Ez_plot);
set(gcf,'Renderer','opengl');
set(h,'AlphaData', ~isnan(Ez_plot));
set(gca,'Color','w');

axis equal tight;
caxis([-2 2]);
title(sprintf('Re{E_z^{sc}}, M = %d', M));
colormap(jet);
colorbar;
