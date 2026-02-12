clc
clear

% Given Parameters
global L a E0 N M beta w d epsilon0 mu0 j

L = 2;
a = 12.5e-3;
f = 300e6;
N = 61;
M = N;

deg90 = pi/2;
deg30 = pi/6;
deg75 = 5*pi/12;

theta_incident = [deg90 deg30 deg75];
theta_plot = linspace(0, 2*pi, N);
E0 = 1;

% Constants
mu0 = 4 * pi * 10^-7;
epsilon0 = 8.854 * 10^-12;
j = 1i;

% variables
w = 2 * pi * f;
beta = w * sqrt(mu0 * epsilon0);
d = L/(N + 1);

% initiate arrays
z_n = zeros(N + 2, 1);
z_m = zeros(M + 2, 1);
v_m = zeros(N, 1);

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

for n_theta = 1:length(theta_incident)

    theta = theta_incident(n_theta);

    % V matrix
    for n = 1:N

        v_m(n) = E0 * exp(j * beta * z_m(n) * cos(theta)) * sin(theta);

        % E_inc,z = -E_theta(0) * sin(theta) * exp(j * beta * z' *
        % cos(theta))
        % v_m = E_inc,v (z'_m)
    
    end

    % I matrix [Z_mn][I_n] = [V_m]
    I_N = Z \ v_m;

    % Radiation
    for n = 1:N

        SRad(n) = s_rad(I_N, theta_plot(n));

    end

% Current Magnitude
figure;
plot(z_n(2:length(z_n) - 1), abs(I_N)*10^3);
xlabel("Z position [m]");
ylabel("Current [mA]");
title("(|I(z)|)")
subtitle("\theta_{inc} = " + theta * 180/pi + "^{\circ}");
grid on;

% Current Phase
figure;
plot(z_n(2:length(z_n) - 1), angle(I_N) * 180 / pi);
xlabel("Z position [m]");
ylabel("Phase [degree]");
title("Phase of I(z)");
subtitle("\theta_{inc} = " + theta * 180/pi + "^{\circ}")
grid on;

fprintf('Maximum |I|_{max} = %.3f mA for \\theta_{inc} = %.3f^{\\circ} and N = %0.3f\n', max(abs(I_N)) * 1e3, theta * 180/pi, N);


end