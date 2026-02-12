% ECE 541 Fall 2025 Project 6
% problem 2 - TM Case

clc
clear

lambda = 1;
a      = 1*lambda;
eps_r  = 4.0;
mu_r   = 1.0;

k  = 2*pi/lambda;
kd = k*sqrt(eps_r*mu_r);

M = ceil(k*a + 15);      

n  = (-M:M).';
jneg = (1j).^(-n);

ka  = k*a;
kda = kd*a;

Jn_ka  = besselj(n,  ka);
Jn_kda = besselj(n,  kda);

Jp_ka  = 0.5*(besselj(n-1,ka)  - besselj(n+1,ka));
Jp_kda = 0.5*(besselj(n-1,kda) - besselj(n+1,kda));

Hn_ka  = besselh(n,2,ka);
Hp_ka  = 0.5*(besselh(n-1,2,ka) - besselh(n+1,2,ka));

den = sqrt(mu_r).*Hp_ka.*Jn_kda - sqrt(eps_r).*Hn_ka.*Jp_kda;

an = -jneg .* ( sqrt(mu_r).*Jp_ka.*Jn_kda - sqrt(eps_r).*Jn_ka.*Jp_kda ) ./ den;

phi = linspace(0, 2*pi, 721);      
ang = phi + pi/2;

E = exp(1j*(n*ang));                 
F = (an.' * E);                       

sigma2D = (4/k) * abs(F).^2;         
sigma2D_norm = sigma2D / lambda;
sigma2D_dB = 10*log10(sigma2D_norm);

figure('Color','w');
plot(rad2deg(phi), sigma2D_dB, 'LineWidth', 1.6);
grid on;
xlim([0 360]);
xlabel('\phi (degrees)');
ylabel('\sigma_{2D}/\lambda (dB)');
title('TM dielectric cylinder: bistatic scattering width');
