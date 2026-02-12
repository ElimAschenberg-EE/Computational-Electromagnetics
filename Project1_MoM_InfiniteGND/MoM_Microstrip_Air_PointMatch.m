% ECE 541 – Applied Electromagnetics – Project 1 (Point Matching)
% w: strip width, h: height, N: segments
% Elim Aschenberg
% 9/27/2025

% Methods of moments: Point Matching

% Define Constants
epsilon0 = 8.8541878188*10^-12;
width = 1;
N = 500;    % # of segments
h = 1;
eta0 = 377;
c0 = 3E8;
epsilon_r = 1;
v1=1;
v2 = 0;

% Initalize variables
delta = width/N;
d = 2*h;

% intialize v and z matrices
v_matrx = v1*ones(1,N);
z = zeros(N,N);

% Populate the matrix z with some values 
for i = 1:N
    xi = delta*i;
    for j = 1:N
        xj = delta*j;
        z(i,j) = 1/(2*pi*epsilon0)*(K1_function(delta/2, xi-xj, 0)-K1_function(delta/2, xi-xj,d));
    end
end

% create charge distribution matrix
Rho = v_matrx * inv(z);

% Plotting
figure;
plot(Rho);
title('Charge Density Distribution for point matching');
xlabel('Index');
ylabel('Rho Value');

% Calculate C'
Charge_Pul = sum(Rho * delta);
Cap_Pul = Charge_Pul/(v1-v2);

% Begin empirical comparison by initalizing matricies
h_val = [0.5,1,2,5,10];
results_vector = zeros(2,length(h_val));
z0_emp = zeros(1, length(h_val));
Cap_Pul_emp = zeros(1, length(z0_emp));
Percent_Error = zeros(2, length (h_val));

for k = 1: length(h_val)
    h = h_val(k);
    ratio = width/h;

        % Empirical Calculations
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
            Cap_Pul_emp (k) = 1/(z0_emp (k) * c0);
        end

% Redefine d for new h values
d = 2*h;

for i = 1:N
    xi = delta*i;
    for j = 1:N
        xj = delta*j;

        % recalculate z for every new h
        z(i,j) = 1/(2*pi*epsilon0)*(K1_function(delta/2, xi-xj, 0)-K1_function(delta/2, xi-xj,d));
    end
end

Rho = v_matrx * inv(z);

Charge_Pul = 0;
for j = 1:N
    Charge_Pul = Charge_Pul + (Rho(j) * delta) ;
end

    % Calculate %error
    Cap_Pul = Charge_Pul/(v1-v2);
    results_vector(1, k) = ratio;
    results_vector(2, k) = Cap_Pul;
    Percent_Error (1, k) = ratio;
    Percent_Error(2, k) = abs(results_vector(2,k)-Cap_Pul_emp(k)) / Cap_Pul_emp(k) * 100;
end

% ----- Plot C' vs w/h -----
figure;
plot(results_vector(1,:), results_vector(2,:),'o-','LineWidth',1.5);
grid on;
xlabel('w/h Ratio');
ylabel('C'' (F/m)');
title(sprintf('Capacitance per Unit Length vs w/h Ratio (N = %d)', N));

% ----- Plot % Error vs w/h -----
figure;
plot(Percent_Error(1,:), Percent_Error(2,:),'o-','LineWidth',1.5);
grid on;
xlabel('w/h Ratio');
ylabel('Relative Percent Error (%)');
title(sprintf('Relative Percentage Error vs w/h Ratio (N = %d)', N));


% ----- Labeled tables -----
fprintf('Capacitance per Unit Length (MoM), N = %d\n', N);
T_Cpul = table(results_vector(1,:).', results_vector(2,:).', ...
    'VariableNames', {'w_over_h','C_Pul'});
disp(T_Cpul);

fprintf('Percent Error, N = %d\n', N);
T_err = table(Percent_Error(1,:).', Percent_Error(2,:).', ...
    'VariableNames', {'w_over_h','PercentError'});
disp(T_err);

% ===== Convergence vs Δ (added) =====
w_fix = width; 
h_fix = h; 
ratio_fix = w_fix/h_fix;

% empirical C' for fixed geometry
if ratio_fix < 1
    p_fix = 0.04*(1-ratio_fix)^2;
else
    p_fix = 0;
end
eps_reff_fix = (epsilon_r + 1)/2 + ((epsilon_r - 1)/2) * ((1 + 12/ratio_fix)^(-1/2) + p_fix);
if ratio_fix > 1
    z0_emp_fix = eta0/sqrt(eps_reff_fix) * (ratio_fix + 1.393 + 0.667*log(ratio_fix + 1.444))^-1;
    C_emp_fix  = 1/(z0_emp_fix * c0);
else
    z0_emp_fix = eta0/(2*pi*sqrt(eps_reff_fix)) * log(8/ratio_fix + ratio_fix/4);
    C_emp_fix  = sqrt(eps_reff_fix)/(z0_emp_fix * c0);
end

% N sweep for Δ = w/N
N_list = round(20 * (1.25).^(0:11)); 
N_list = unique(N_list);

% compute C', error for each Δ
Delta_pm = zeros(size(N_list));
C_pm     = zeros(size(N_list));
Err_pm   = zeros(size(N_list));
for t = 1:numel(N_list)
    N_try = N_list(t);
    Delta_pm(t) = w_fix / N_try;
    C_pm(t) = local_Cprime_pointmatch(N_try, w_fix, h_fix, epsilon0);
    Err_pm(t) = abs(C_pm(t) - C_emp_fix);
end

% sort errors (desc) and compute rates
[Err_sorted_pm, idx_pm] = sort(Err_pm, 'descend');
Delta_sorted_pm = Delta_pm(idx_pm);
N_sorted_pm     = N_list(idx_pm);
C_sorted_pm     = C_pm(idx_pm);
Rate_PM = NaN(size(Err_sorted_pm));
for i = 1:numel(Err_sorted_pm)-1
    Rate_PM(i) = log2(Err_sorted_pm(i)/Err_sorted_pm(i+1)) / log2(Delta_sorted_pm(i)/Delta_sorted_pm(i+1));
end

% table + log–log plot
T_PM = table(N_sorted_pm.', Delta_sorted_pm.', C_sorted_pm.', Err_sorted_pm.', Rate_PM.', ...
    'VariableNames', {'N','Delta','C_Pul_PM','AbsError','Rate'});
disp(T_PM);

figure;
loglog(Delta_pm, Err_pm, 'o-','LineWidth',1.5); grid on;
xlabel('\Delta = w/N');
ylabel('|C''_{PM} - C''_{emp}|');
title('Point Matching Convergence: Error vs \Delta');

% ===== helper for convergence (added) =====
function Cprime = local_Cprime_pointmatch(N, w, h, epsilon0)
    % build Z, solve rho, return C' for given N
    v1 = 1; v2 = 0;
    delta = w/N; 
    d = 2*h;
    zloc = zeros(N,N);
    for i = 1:N
        xi = delta*i;
        for j = 1:N
            xj = delta*j;
            zloc(i,j) = 1/(2*pi*epsilon0)*(K1_function(delta/2, xi-xj, 0)-K1_function(delta/2, xi-xj, d));
        end
    end
    Rho_loc = (ones(1,N))/zloc;
    Q = sum(Rho_loc)*delta;
    Cprime = Q/(v1 - v2);
end
