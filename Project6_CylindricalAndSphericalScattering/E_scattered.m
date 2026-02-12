% ECE 541 Fall 2025 Project 6
% problem 2 - Scattered E

clc
clear

lambda = 1;                
a      = 1*lambda;         
eps_r  = 4.0;               
mu_r   = 1.0;               

k  = 2*pi/lambda; 
kd = k*sqrt(eps_r*mu_r);

M = 40;

x = -3*lambda : 0.02*lambda : 3*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);

Einc = exp(-1j*k*X);  


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

cn = (1j).^(-(n+1)) .* (2*sqrt(mu_r)) ./ (pi*ka .* den);

Esc_out = zeros(size(RHO));
Eint    = zeros(size(RHO));

for ii = 1:numel(n)
    ni = n(ii);
    phase = exp(1j*ni*PHI);

    Esc_out = Esc_out + an(ii) .* besselh(ni,2,k*RHO)  .* phase;

    Eint    = Eint    + cn(ii) .* besselj(ni,  kd*RHO) .* phase;
end

Esc = Esc_out;
inside = (RHO < a);
Esc(inside) = Eint(inside) - Einc(inside);

% ---------------- Plot ----------------
figure('Color','w');
imagesc(x/lambda, y/lambda, real(Esc)); 
axis equal tight;
set(gca,'YDir','normal');
xlabel('x/\lambda'); ylabel('y/\lambda');
caxis([-2 2]);
title('TM dielectric cylinder: Re\{E_z^{sc}/E_0\}');
colormap(jet)
colorbar;

% Cylinder boundary
hold on;
th = linspace(0,2*pi,600);
plot((a/lambda)*cos(th),(a/lambda)*sin(th),'k','LineWidth',1.2);
hold off;
