% ECE 541 Fall 2025 Project 6
% problem 1(c) - TE scattering width

clc
clear

j = 1i;
lambda = 1;
a = 1*lambda;
k = 2*pi/lambda;
ka = k*a;

M = 40;
n = (-M:M).';      

c_n = j.^(-n);

J_ka_prime = 0.5*( besselj(n-1,ka) - besselj(n+1,ka) );
H_ka_prime = 0.5*( besselh(n-1,2,ka) - besselh(n+1,2,ka) );

b_n = -c_n .* ( J_ka_prime ./ H_ka_prime );

phi = linspace(0, 2*pi, 2000);
F = ( (b_n .* (j.^n)).' * exp(1i*n*phi) );

sigma = (4/k) * abs(F).^2;

Sigma_total = trapz(phi, sigma);

% Plot
figure;
plot(phi*180/pi, sigma, 'LineWidth', 1.2);
grid on;
xlim([0 360]);
xticks(0:30:360);
xlabel('\phi (deg)');
ylabel('\sigma(\phi)  (length)');
title(sprintf('TE_z PEC Cylinder Scattering Width, ka=%.3f, M=%d', ka, M));
[~,i0]   = min(abs(phi - 0));
[~,i180] = min(abs(phi - pi));
fprintf('Forward scatter (phi=0 deg):   sigma_f = %.6g\n', sigma(i0));
fprintf('Backscatter (phi=180 deg):     sigma_b = %.6g\n', sigma(i180));
fprintf('Total scattering width:        Sigma   = %.6g\n', Sigma_total);
