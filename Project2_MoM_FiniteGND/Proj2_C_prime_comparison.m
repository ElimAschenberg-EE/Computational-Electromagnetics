clc; clear;

eps0 = 8.854187817e-12;
er   = 1;
eps  = eps0 * er;
coef = 1/(2*pi*eps);
v1 = 1;
v2 = 0;

w1_list = [0.5, 1.0, 2.0];      
h_list  = [0.5, 1.0, 2.0];      
N1_list = [50, 100, 200];       % Number of top subdivisions
% w2 will automatically be set to w1 + 10*h (per project requirement)

max_rows = numel(w1_list)*numel(h_list)*numel(N1_list);
results  = zeros(max_rows, 9);
row = 1;

% ============================================================
%                      MAIN CALCULATION
% ============================================================
for ia = 1:numel(w1_list)
  for ic = 1:numel(h_list)

    % ----- Geometry -----
    w1 = w1_list(ia);
    h  = h_list(ic);
    w2_target = w1 + 10*h;   % Large ground width per requirement

    for kn = 1:numel(N1_list)

      % ----- Discretization -----
      N1 = N1_list(kn);
      delta = w1 / N1;
      a = delta / 2;
      N2 = max(1, round(w2_target / delta));
      w2 = N2 * delta;
      Ntot = N1 + N2;

      % ============================================================
      %        APPROACH I — POINT MATCHING (Pulse–Point)
      % ============================================================

      % --- Build Kernels (fast via offset indexing) ---
      ker_tt_I = zeros(1,N1);
      for k = 1:N1
        dx = (k-1)*delta;
        ker_tt_I(k) = coef * K1_function(a, dx, 0);
      end
      ker_gg_I = zeros(1,N2);
      for k = 1:N2
        dx = (k-1)*delta;
        ker_gg_I(k) = coef * K1_function(a, dx, 0);
      end
      L = N1 + N2 - 1;
      ker_tg_I = zeros(1,L);
      for k = 1:L
        dx = (k - N2)*delta;
        ker_tg_I(k) = coef * K1_function(a, dx, h);
      end

      % --- Build Impedance Matrix Z_I ---
      Ztt_I = toeplitz(ker_tt_I, ker_tt_I);
      Zgg_I = toeplitz(ker_gg_I, ker_gg_I);
      Ztg_I = zeros(N1,N2);
      for i = 1:N1
        for j = 1:N2
          o = i - j;
          Ztg_I(i,j) = ker_tg_I(o + N2);
        end
      end
      Zgt_I = Ztg_I.';
      Z_I = [Ztt_I, Ztg_I; Zgt_I, Zgg_I];

      % --- Gauge Fixing and Solve ---
      Z_I = Z_I - Z_I(1,:);
      Z_I(1,:) = delta;
      V_I = [zeros(N1,1); (v2 - v1)*ones(N2,1)];
      rho_I = Z_I \ V_I;

      Qp_I = sum(rho_I(1:N1)) * delta;
      Cprime_I = Qp_I / (v1 - v2);

      % ============================================================
      %        APPROACH II — GALERKIN (Pulse–Pulse)
      % ============================================================

      % --- Build Kernels (fast) ---
      ker_tt_II = zeros(1,N1);
      for k = 1:N1
        dx = (k-1)*delta;
        ker_tt_II(k) = coef * K2_function(a, dx + a, 0, dx - a, 0);
      end
      ker_gg_II = zeros(1,N2);
      for k = 1:N2
        dx = (k-1)*delta;
        ker_gg_II(k) = coef * K2_function(a, dx + a, 0, dx - a, 0);
      end
      ker_tg_II = zeros(1,L);
      for k = 1:L
        dx = (k - N2)*delta;
        ker_tg_II(k) = coef * K2_function(a, dx + a, h, dx - a, h);
      end

      % --- Build Impedance Matrix Z_II ---
      Ztt_II = toeplitz(ker_tt_II, ker_tt_II);
      Zgg_II = toeplitz(ker_gg_II, ker_gg_II);
      Ztg_II = zeros(N1,N2);
      for i = 1:N1
        for j = 1:N2
          o = i - j;
          Ztg_II(i,j) = ker_tg_II(o + N2);
        end
      end
      Zgt_II = Ztg_II.';
      Z_II = [Ztt_II, Ztg_II; Zgt_II, Zgg_II];

      % --- Gauge Fixing and Solve ---
      Z_II = Z_II - Z_II(1,:);
      Z_II(1,:) = delta;
      V_II = delta * [zeros(N1,1); (v1 - v2)*ones(N2,1)];
      rho_II = Z_II \ V_II;

      Qp_II = sum(rho_II(1:N1)) * delta;
      Cprime_II = Qp_II / (v1 - v2);

      % ============================================================
      %        COMPARISON AND STORAGE
      % ============================================================
      RelErr = 100 * abs(Cprime_I - Cprime_II) / max(eps, abs(Cprime_II));

      results(row,:) = [w1, w2, h, N1, N2, Ntot, Cprime_I, Cprime_II, RelErr];
      row = row + 1;

    end
  end
end

results = results(1:row-1,:);

% ============================================================
%                        TABULAR OUTPUT
% ============================================================
T = table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), results(:,6), ...
          results(:,7), results(:,8), results(:,9), ...
   'VariableNames', {'w1_m','w2_used_m','h_m','N1','N2','N_total', ...
                     'Cprime_PointMatching_Fperm','Cprime_Galerkin_Fperm','RelErr_percent'});

disp(' ');
disp('============================================================');
disp('  COMPARISON OF POINT MATCHING AND GALERKIN APPROACHES  ');
disp('============================================================');
disp(T);

% Subset where w2 - w1 > 10*h (should all satisfy but we'll confirm)
mask_far = (results(:,2) - results(:,1)) > (10*results(:,3));
T_far = T(mask_far,:);
disp(' ');
disp('============================================================');
disp('  CASES WHERE w2 - w1 > 10h  (FAR-GROUND CONDITION)');
disp('============================================================');
disp(T_far);

% ============================================================
%                        PROFESSIONAL PLOTS
% ============================================================
% Plot 1: Comparison of C' for both methods vs N1
triplets = unique(results(:,[1,3]), 'rows', 'stable');  % unique (w1,h)
figure; hold on; grid on;
for t = 1:size(triplets,1)
    w1c = triplets(t,1);
    hc  = triplets(t,2);
    mask = (results(:,1)==w1c) & (results(:,3)==hc);
    sub = results(mask,:);
    [Nvec, idx] = sort(sub(:,4));
    C1 = sub(idx,7);  % Point Matching
    C2 = sub(idx,8);  % Galerkin
    plot(Nvec, C1*1e12, '-o', 'LineWidth', 1.2, 'DisplayName', ...
        sprintf('Point Matching: w1=%.2g, h=%.2g', w1c, hc));
    plot(Nvec, C2*1e12, '--s', 'LineWidth', 1.2, 'DisplayName', ...
        sprintf('Galerkin: w1=%.2g, h=%.2g', w1c, hc));
end
xlabel('Number of Top Conductor Segments, N_1');
ylabel('Capacitance per Unit Length, C''  (pF/m)');
title('Comparison of Capacitance per Unit Length (Point Matching vs Galerkin)');
legend('Location','northeastoutside');

% Plot 2: Relative % Error between the Two Approaches
figure; hold on; grid on;
for t = 1:size(triplets,1)
    w1c = triplets(t,1);
    hc  = triplets(t,2);
    mask = (results(:,1)==w1c) & (results(:,3)==hc);
    sub = results(mask,:);
    [Nvec, idx] = sort(sub(:,4));
    E = sub(idx,9);
    plot(Nvec, E, '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('w1=%.2g, h=%.2g', w1c, hc));
end
xlabel('Number of Top Conductor Segments, N_1');
ylabel('Relative Percent Error, |C''_{PM} - C''_{Gal}| / C''_{Gal} × 100');
title('Relative Percentage Error Between Point Matching and Galerkin Methods');
legend('Location','northeastoutside');
