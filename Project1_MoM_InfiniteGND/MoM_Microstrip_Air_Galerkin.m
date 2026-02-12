% ECE 541 – Applied Electromagnetics – Project 1 (Galerkin)
% w: strip width, h: height, N: segments
% Elim Aschenberg
% 9/27/2025

% Galerkin Method

% initalizing constants
epsilon0 = 8.854187817e-12;
w = 1; 
h = 1; 
N = 1000;
v1 = 1;
v2 = 0;
epsilon_r = 1;
eta0 = 377;
c0 = 3E8;

% initalizing variables
d = 2 * h;
delta = w/N;
a = delta/2;

% creating matricies
x = zeros(N,1);
v = delta*ones(N,1);
z = zeros(N,N);

% calculate x positions at the center of each segment
for i = 1:N
    x(i) = (i-0.5) * delta;
end

for i = 1:N
    xi = x(i);
    for j = 1:N
        xj = x(j);
        x_right = xi - (xj + a);
        x_left  = xi - (xj - a);
        K_self  = K2_function(a, x_right, 0, x_left, 0);
        K_image = K2_function(a, x_right, d, x_left, d);
        z(i,j) = -1/(2*pi*epsilon0) * (K_self - K_image);
    end
end

% Z(i,j) = Z(j,i) -> matrix is symmetric
for i = 1:N
    for j = 1:N
        z(i,j) = 0.5 * (z(i,j) + z(j,i));
    end
end

% Calculations for C'
Rho = z \ v;
Charge_Pul = sum(Rho * delta);
Cap_Pul = Charge_Pul/(v1-v2);

% Plotting
figure; 
plot(x, Rho/max(abs(Rho)),'-'); 
grid on
xlabel('x'); 
ylabel('\rho(x)'); 
title('Charge Density Distribution');
drawnow;

% Setting up matricies
h_val = [0.5,1,2,5,10];
results_vector = zeros(2,length(h_val));
z0_emp = zeros(1, length(h_val));
Cap_Pul_emp = zeros(1, length(z0_emp));
Percent_Error = zeros(2, length (h_val));

for k = 1: length(h_val)
    h = h_val(k);

    % Redefine variables reliant upon h
    ratio = w/h;
    d = 2 * h;

    % empirical calulations
    if ratio < 1
        p = 0.04*(1-ratio)^2;
    else
        p = 0;
    end
    epsilon_reff = (epsilon_r + 1)/2 + ((epsilon_r - 1)/2) * ((1 + 12/ratio)^-(1/2) + p);
    if ratio > 1
        z0_emp (k) = eta0/sqrt(epsilon_reff) * (ratio + 1.393 + 0.667 * log (ratio + 1.444))^-1;
        Cap_Pul_emp (k) = 1/(z0_emp (k) * c0);
    else
        z0_emp (k) = eta0/(2*pi*sqrt(epsilon_reff)) * log (8/ratio + ratio/4);
        Cap_Pul_emp (k) = sqrt(epsilon_reff)/(z0_emp (k) * c0);
    end

    % Repeating the Galerkin code above
    for i = 1:N
        x(i) = (i-0.5) * delta;
    end
    for i = 1:N
        xi = x(i);
        for j = 1:N
            xj = x(j);
            x_right = xi - (xj + a);
            x_left  = xi - (xj - a);
            K_self  = K2_function(a, x_right, 0, x_left, 0);
            K_image = K2_function(a, x_right, d, x_left, d);
            z(i,j) = -1/(2*pi*epsilon0) * (K_self - K_image);
        end
    end
    for i = 1:N
        for j = 1:N
            z(i,j) = 0.5 * (z(i,j) + z(j,i));
        end
    end

    Rho = z \ v;

    Charge_Pul = 0;
    for j = 1:N
        Charge_Pul = Charge_Pul + (Rho(j) * delta);
    end

    Cap_Pul = Charge_Pul/(v1-v2);
    results_vector(1, k) = ratio;
    results_vector(2, k) = Cap_Pul;
    Percent_Error (1, k) = ratio;
    Percent_Error(2, k) = abs(results_vector(2,k)-Cap_Pul_emp(k))...
                            / Cap_Pul_emp(k) * 100;
end

% %Error plot with N in title
figure;
plot(Percent_Error(1,:), Percent_Error(2,:),'o-','LineWidth',1.5); grid on;
xlabel('w/h Ratio');
ylabel('Relative Percent Error (%)');
title(sprintf('Relative Percentage Error vs w/h Ratio (N = %d)', N));
drawnow;

% C' vs w/h with N in title
figure;
plot(results_vector(1,:), results_vector(2,:),'o-','LineWidth',1.5); grid on;
xlabel('w/h Ratio');
ylabel('C'' (F/m)');
title(sprintf('Capacitance per Unit Length vs w/h Ratio (N = %d)', N));
drawnow;

% labeled tables with N in header
fprintf('Capacitance per Unit Length (MoM), N = %d\n', N);
T_Cpul = table(results_vector(1,:).', results_vector(2,:).', ...
    'VariableNames', {'w_over_h','C_Pul'});
disp(T_Cpul);

fprintf('Percent Error, N = %d\n', N);
T_err = table(Percent_Error(1,:).', Percent_Error(2,:).', ...
    'VariableNames', {'w_over_h','PercentError'});
disp(T_err);

% Convergence 
if ratio < 1
        p = 0.04*(1-ratio)^2;
    else
        p = 0;
end

    epsilon_reff = (epsilon_r + 1)/2 + ((epsilon_r - 1)/2) * ((1 + 12/ratio)^-(1/2) + p);

    if ratio > 1
        z0_emp (k) = eta0/sqrt(epsilon_reff) * (ratio + 1.393 + 0.667 * log (ratio + 1.444))^-1;
        Cap_Pul_emp (k) = 1/(z0_emp (k) * c0);
    else
        z0_emp (k) = eta0/(2*pi*sqrt(epsilon_reff)) * log (8/ratio + ratio/4);
        Cap_Pul_emp (k) = sqrt(epsilon_reff)/(z0_emp (k) * c0);
    end

N_list = 20:20:260;   % multiples of 20 up to 260

delta_conv = zeros(size(N_list));
C_conv = zeros(size(N_list));
Err_abs = zeros(size(N_list));

for t = 1:numel(N_list)
    N_try = N_list(t);
    delta_conv(t) = w / N_try;
    C_conv(t) = local_Cprime_galerkin(N_try, w, h, epsilon0, epsilon_r);
    Err_abs(t) = abs(C_conv(t) - C_emp_fix);
end

% sort errors (desc) and compute rates
[Err_sorted, idx] = sort(Err_abs, 'descend');
Delta_sorted = delta_conv(idx);
N_sorted     = N_list(idx);
C_sorted     = C_conv(idx);
Rate = NaN(size(Err_sorted));
for i = 1:numel(Err_sorted)-1
    Rate(i) = log2(Err_sorted(i)/Err_sorted(i+1)) / log2(Delta_sorted(i)/Delta_sorted(i+1));
end

% table + log–log plot
T_GAL = table(N_sorted.', Delta_sorted.', C_sorted.', Err_sorted.', Rate.', ...
    'VariableNames', {'N','Delta','C_Pul_Galerkin','AbsError','Rate'});
disp(T_GAL);

figure;
loglog(delta_conv, Err_abs, 'o-','LineWidth',1.5); grid on;
xlabel('\Delta = w/N');
ylabel('|C''_{GAL} - C''_{emp}|');
title('Galerkin Convergence: Error vs \Delta');

% ===== helper for convergence (added) =====
function Cprime = local_Cprime_galerkin(N, w, h, epsilon0, epsilon_r)
    % build Z, solve rho, return C' for given N
    v1 = 1; 
    v2 = 0; 
    eta0 = 377; 
    c0 = 3e8;
    d = 2*h; 
    delta = w/N; 
    a = delta/2;

    % centers
    x = ((1:N)-0.5).' * delta;

    % Z matrix
    zmat = zeros(N,N);
    for i = 1:N
        xi = x(i);
        for j = 1:N
            xj = x(j);
            x_right = xi - (xj + a);
            x_left  = xi - (xj - a);
            K_self  = K2_function(a, x_right, 0, x_left, 0);
            K_image = K2_function(a, x_right, d, x_left, d);
            zmat(i,j) = -1/(2*pi*epsilon0) * (K_self - K_image);
        end
    end
    zmat = 0.5*(zmat + zmat.');

    % solve and C'
    vvec = delta*ones(N,1);
    rho  = zmat \ vvec;
    Q    = sum(rho) * delta;
    Cprime = Q/(v1 - v2);
end
