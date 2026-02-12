clc
clear

% Given Parameters
global L a E0 N M beta w d epsilon0 mu0 j

f = 300e6;
N = 50;
M = N;
L = linspace(0.25, 2, N);

theta = pi/2;
theta_plot = linspace(0, 2*pi, N);
E0 = 1;

% Constants
mu0 = 4 * pi * 10^-7;
epsilon0 = 8.854 * 10^-12;
j = 1i;

% variables
w = 2 * pi * f;
beta = w * sqrt(mu0 * epsilon0);

% initiate arrays
z_n = zeros(N + 2, 1);
z_m = zeros(M + 2, 1);
v_m = zeros(N, 1);

for i = 1:N

    d = L(i)/(N + 1);
    a = L(i)/160;

    % fill z and z' matrices with locations
    for n = 1:N 

        z_n(n + 1) = z_n(n) + d;
        z_m(n + 1) = z_m(n) + d;

    end

    for m = 1:M

        for n = 1:N

            Z(m,n) = j * w * mu0 / (4 * pi) * (psi_funct(1, n, 0, m) + 1 / (beta^2 * d) * psi_funct(2, n, 0, m));
        
        end

    end


    % V matrix
    for n = 1:N

        v_m(n) = -E0 * exp(j * beta * z_m(n) * cos(theta)) * sin(theta);

        % E_inc,z = -E_theta(0) * sin(theta) * exp(j * beta * z' *
        % cos(theta))
        % v_m = E_inc,v (z'_m)
    
    end

    % I matrix [Z_mn][I_n] = [V_m]
    I_N = Z \ v_m;



    SRad(i) = s_rad(I_N, pi/2);



    % Plot Radiation pattern cross section
    %{
    figure;
    polarplot(abs(SRad));
    hold on;
    polarplot([theta; theta], [0; max(abs(SRad))]);
    legend("S_{rad}", "\Theta_{inc}");
    title("Length = " + L(i) + ", radius =" + a + " & \theta_{inc} = " + theta * 180/pi + " degrees");
    hold off;
    %}

end

figure;
plot(L, SRad, '-o');
xlabel('Length (m)');
ylabel('radiation magnitude at \theta_{inc} = 90^{\circ} (U)')
title("Monostatic radar cross-section of the scatterer for L ranging from 0.25m to 2m");

display(SRad)
display(L)
Display(a)