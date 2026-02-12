% Colorado State University
% ECE 541 – Applied Electromagnetics – Project 2
% Elim Aschenberg
% 10/9/2025

% MoM Point Matching Finite Ground Plane

clc; clear;

eps0 = 8.8541878188e-12;
er   = 1;
eps  = eps0 * er;
coef = 1/(2*pi*eps);
eta0 = 377;
c0   = 3e8;
v1 = 1;
v2 = 0;
w1 = 1; 
w2 = 20; 
N1 = 20; 
h = 1;

delta = w1 / N1;
a = delta/2;
N2 = max(1, round(w2 / delta));
w2 = N2 * delta;
N  = N1 + N2;

x1 = linspace(-w1/2 + a, +w1/2 - a, N1).';  y1 = h*ones(N1,1);
x2 = linspace(-w2/2 + a, +w2/2 - a, N2).';  y2 = zeros(N2,1);

Z1_self = zeros(N1,N1);
for i = 1:N1
  for j = 1:N1
    Z1_self(i,j) = coef * K1_function(a, x1(i)-x1(j), y1(i)-y1(j));
  end
end

Z2_self = zeros(N2,N2);
for i = 1:N2
  for j = 1:N2
    Z2_self(i,j) = coef * K1_function(a, x2(i)-x2(j), y2(i)-y2(j));
  end
end

Z1_to_Z2 = zeros(N1,N2);
for i = 1:N1
  for j = 1:N2
    Z1_to_Z2(i,j) = coef * K1_function(a, x1(i)-x2(j), y1(i)-y2(j));
  end
end

Z2_to_Z1 = zeros(N2,N1);
for i = 1:N2
  for j = 1:N1
    Z2_to_Z1(i,j) = coef * K1_function(a, x2(i)-x1(j), y2(i)-y1(j));
  end
end

% Make Z matrix
Z = [Z1_self, Z1_to_Z2; Z2_to_Z1, Z2_self];
Z = Z - Z(1,:);
Z(1,:) = delta;

% Make V matrix
V = [zeros(N1,1); (v2 - v1)*ones(N2,1)];

rho = Z \ V;
rho1 = rho(1:N1);
rho2 = rho(N1+1:end);

Q_prime = sum(rho1)*delta;
C_prime = Q_prime/(v1 - v2);

% Plot charge distribution
fig_dist = figure; hold on; grid on;
plot(x1, rho1, 'LineWidth', 1.2);
plot(x2, rho2, 'LineWidth', 1.2);
xlabel('x position (m)'); ylabel('\rho_s (C/m^2)');
legend('Conductor (y=h)','Ground (y=0)','Location','best');
title(sprintf('Finite Ground Charge Distribution  w_1=%.1f, w_2=%.1f, h=%.1f', w1, w2, h));

% Create multiple sets of w1,w2,h, and N
w1_list = [0.5, 1.0];
w2_list = [2.0, 10.0];
h_list  = [0.5, 1.0, 2.0];
N1_list = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]; % your custom N1s

% Figure 1: % Error vs N1
fig_err  = figure; hold on; grid on;
xlabel('N_1 (Segments on conductor)');
ylabel('Relative %Error vs Empirical');
title('Relative %Error (Empirical) vs N_1');

% Figure 2: Convergence vs N1
fig_conv = figure; hold on; grid on;
xlabel('N_1 (Segments on conductor)');
ylabel('Relative Error to C''_{ref} (%)');
title('Convergence to C''_{ref} vs N_1 (log scale)');

% Results table to help with reporting
max_rows = length(w1_list)*length(w2_list)*length(h_list)*length(N1_list);
results = zeros(max_rows, 9);
row = 1;

for ia = 1:length(w1_list)
  for ib = 1:length(w2_list)
    for ic = 1:length(h_list)

      w1 = w1_list(ia);
      w2_target = w2_list(ib);
      h = h_list(ic);

      % Empirical C' calculations
      ratio = w1 / h;
      if ratio < 1, p = 0.04*(1 - ratio)^2; else, p = 0; end
      eps_eff = (er + 1)/2 + (er - 1)/2 * ((1 + 12/ratio)^(-0.5) + p);
      if ratio <= 1
        Z0_emp = (eta0/(2*pi*sqrt(eps_eff))) * log(8/ratio + ratio/4);
      else
        Z0_emp = (eta0/sqrt(eps_eff)) / (ratio + 1.393 + 0.667*log(ratio + 1.444));
      end
      C_emp = 1/(Z0_emp * c0);

      % Arrays for this geometry
      Percent_error_empirical = zeros(size(N1_list));
      Cap_point_matching_list = zeros(size(N1_list));

      % Sweep N1
      for k = 1:length(N1_list)
        N1 = N1_list(k);
        delta = w1 / N1;
        a = delta/2;

        N2 = max(1, round(w2_target / delta));
        w2_current = N2 * delta;

        % Centers with new geometry
        x1 = linspace(-w1/2 + a, +w1/2 - a, N1).';  yt = h*ones(N1,1);
        x2 = linspace(-w2_current/2 + a, +w2_current/2 - a, N2).';  yg = zeros(N2,1);

        Z1_self = zeros(N1,N1);
        for i = 1:N1
          for j = 1:N1
            Z1_self(i,j) = coef * K1_function(a, x1(i)-x1(j), yt(i)-yt(j));
          end
        end

        Z2_self = zeros(N2,N2);
        for i = 1:N2
          for j = 1:N2
            Z2_self(i,j) = coef * K1_function(a, x2(i)-x2(j), yg(i)-yg(j));
          end
        end

        Z1_to_Z2 = zeros(N1,N2);
        for i = 1:N1
          for j = 1:N2
            Z1_to_Z2(i,j) = coef * K1_function(a, x1(i)-x2(j), yt(i)-yg(j));
          end
        end

        Z2_to_Z1 = zeros(N2,N1);
        for i = 1:N2
          for j = 1:N1
            Z2_to_Z1(i,j) = coef * K1_function(a, x2(i)-x1(j), yg(i)-yt(j));
          end
        end

        % Make new Z matrix
        Z = [Z1_self, Z1_to_Z2; Z2_to_Z1, Z2_self];
        Z = Z - Z(1,:);
        Z(1,:) = delta;

        % Re-implement V matrix
        V = [zeros(N1,1); (v2 - v1)*ones(N2,1)];
        rho = Z \ V;

        % C' and empirical % error
        Q_prime = sum(rho(1:N1)) * delta;
        C_prime = Q_prime / (v1 - v2);
        Cap_point_matching_list(k) = C_prime;
        Percent_error_empirical(k) = 100 * abs(C_prime - C_emp) / C_emp;

        % Save a row for reporting (note: w2_current depends on N1)
        results(row,:) = [w1, w2_current, h, N1, N2, N1+N2, C_prime, C_emp, Percent_error_empirical(k)];
        row = row + 1;
      end

      % Graph 1: %error vs N1
      figure(fig_err.Number);
      plot(N1_list, Percent_error_empirical, '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('w1=%.2g, w2=%.2g, h=%.2g', w1, w2_target, h));

      % Graph 2: Convergence vs N1
      C_ref = Cap_point_matching_list(end);
      % use machine epsilon, not your 'eps' variable
      ErrToRefPct = 100 * abs(Cap_point_matching_list - C_ref) / max(eps(1.0), abs(C_ref));
      figure(fig_conv.Number);
      semilogy(N1_list, ErrToRefPct, '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('w1=%.2g, w2=%.2g, h=%.2g', w1, w2_target, h));

    end
  end
end

% Finish figures
figure(fig_err.Number);
yline(1,'--','1% target','LabelHorizontalAlignment','left');
legend('show','Location','northeastoutside');

figure(fig_conv.Number);
legend('show','Location','northeastoutside');

% Table
results = results(1:row-1, :);
T = table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), results(:,6), ...
          results(:,7), results(:,8), results(:,9), ...
          'VariableNames', {'w1_m','w2_used_m','h_m','N1','N2','N_total', ...
                            'Cprime_MoM_Fperm','Cprime_Emp_Fperm','Error_percent'});
disp(' ');
disp('All runs (each row = one geometry at one N1):');
disp(T);

%% ---------- Point Matching: convergence RATE vs Delta (assignment formula) ----------
% r_i = log2(err_i/err_{i+1}) / log2(Δ_i/Δ_{i+1}),  err_i = |C'(N1_i) - C_inf|,  Δ = w1/N1
% Group by (w1,h) AND by rounded w2_used to avoid mixing different grounds.

triplets_wh = unique(results(:,[1,3]), 'rows', 'stable');  % [w1, h]

figure('Name','Point Matching: convergence rate vs Delta'); hold on; grid on;
xlabel('\Delta = w_1/N_1 (m)'); set(gca,'XScale','log'); set(gca,'XDir','reverse');
ylabel('Convergence rate  r_i');
title('Point Matching convergence rate  r_i = log_2(err_i/err_{i+1}) / log_2(\Delta_i/\Delta_{i+1})');

for t = 1:size(triplets_wh,1)
    w1c = triplets_wh(t,1);  hc = triplets_wh(t,2);

    mask_wh = (results(:,1)==w1c) & (results(:,3)==hc);
    sub_wh  = results(mask_wh,:);  % [w1 w2_used h N1 N2 Ntot C_MoM C_emp Err%]
    if isempty(sub_wh), continue; end

    % sub-group by rounded w2_used so each curve corresponds to one ground width
    w2_bins = unique(round(sub_wh(:,2), 6));   % round to 1e-6 m
    for b = 1:numel(w2_bins)
        w2b = w2_bins(b);
        sub  = sub_wh( round(sub_wh(:,2),6) == w2b, : );
        if size(sub,1) < 3, continue; end

        % sort by N1, compute Δ, C_inf, errors
        [Nvec, idx] = sort(sub(:,4));     % N1 ascending
        Cvec = sub(idx,7);
        Dlt  = w1c ./ Nvec;               % Δ = w1/N1

        % C_inf from two finest meshes
        N1a=Nvec(end-1); N1b=Nvec(end); C1a=Cvec(end-1); C1b=Cvec(end);
        Cinf = (N1b*C1b - N1a*C1a)/(N1b - N1a);

        Err = abs(Cvec - Cinf);

        % Δ descending (coarse→fine)
        [Ddesc, idd] = sort(Dlt,'descend');  Edesc = Err(idd);
        if numel(Ddesc) < 2, continue; end

        % pairwise rates and plotting abscissa (geometric-mean Δ)
        numPairs = numel(Ddesc)-1;
        rates = log2(Edesc(1:numPairs)./Edesc(2:numPairs+1)) ./ ...
                log2(Ddesc(1:numPairs)./Ddesc(2:numPairs+1));
        xpair = sqrt(Ddesc(1:numPairs) .* Ddesc(2:numPairs+1));

        finiteMask = isfinite(rates) & isfinite(xpair);
        if any(finiteMask)
            semilogx(xpair(finiteMask), rates(finiteMask), '-o', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('w1=%.2g, w2\\approx%.3g, h=%.2g', w1c, w2b, hc));

            % (optional) console printout
            fprintf('\n[PM] (w1=%.2g, h=%.2g; w2≈%.6g)\n', w1c, hc, w2b);
            fprintf('  i    N1_i   N1_{i+1}        Δ_i          Δ_{i+1}      r_i\n');
            for i=1:numPairs
                fprintf('%3d  %6d   %8d   %12.5e  %12.5e   %7.3f\n', ...
                    i, Nvec(idd(i)), Nvec(idd(i+1)), Ddesc(i), Ddesc(i+1), rates(i));
            end
        end
    end
end
legend('show','Location','bestoutside'); hold off;
