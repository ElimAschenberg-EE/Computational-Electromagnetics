% ECE 541 Fall 2025 Project 6
% problem 1(c) - TM scattering width

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
b_n = -c_n .* ( besselj(n,ka) ./ besselh(n,2,ka) );

phi = linspace(0,2*pi,2000);
F = ( (b_n .* (j.^n)).' * exp(1i*n*phi) ); 

sigma = (4/k) * abs(F).^2;

Sigma_total = trapz(phi, sigma);

% Plot
figure;
plot(phi*180/pi, sigma, 'LineWidth', 1.2);
xlim([0 360]);
xticks(0:30:360);
grid on;
xlabel('\phi (deg)');
ylabel('\sigma(\phi)  (length)');
title(sprintf('TM_z PEC Cylinder Scattering Width, ka=%.3f, M=%d', ka, M));
fprintf('Backscatter (phi=180 deg): sigma_b = %.6g\n', sigma(round(end/2)));
fprintf('Forward scatter (phi=0 deg): sigma_f = %.6g\n', sigma(1));
fprintf('Total scattering width (integrated): Sigma = %.6g\n', Sigma_total);
