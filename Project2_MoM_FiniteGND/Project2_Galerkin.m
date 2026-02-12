% Colorado State University
% ECE 541 – Applied Electromagnetics – Project 2
% Elim Aschenberg
% 10/9/2025

% MoM Galerkin Finite Ground Plane

clc; clear;

epsilon0 = 8.854187817e-12;
epsilon_r = 1;
epsilon = epsilon0*epsilon_r;
w1 = 1;                 
w2 = 10;                 
h  = 1;                 
v1 = 1; 
v2 = 0;         
N1 = 100;
eta0 = 377;
c0 = 3e8;

delta = w1 / N1;
a = delta/2;
N2 = max(1, round(w2 / delta));
w2 = N2 * delta;
N = N1 + N2;

coef = 1/(2*pi*epsilon);

% set conductor centered above ground 
x1 = linspace(-w1/2 + a, +w1/2 - a, N1).';  y1 = h*ones(N1,1);
x2 = linspace(-w2/2 + a, +w2/2 - a, N2).';  y2 = zeros(N2,1);

Z1_self = zeros(N1,N1);
for i = 1:N1
    for j = 1:N1
        dx = x1(i) - x1(j);
        x_right = dx + a;  
        x_left = dx - a;
        Z1_self(i,j) = coef * K2_function(a, x_right, 0, x_left, 0);
    end
end

Z2_self = zeros(N2,N2);
for i = 1:N2
    for j = 1:N2
        dx = x2(i) - x2(j);
        x_right = dx + a;  
        x_left = dx - a;
        Z2_self(i,j) = coef * K2_function(a, x_right, 0, x_left, 0);
    end
end

Z1_to_Z2 = zeros(N1,N2);
for i = 1:N1
    for j = 1:N2
        dx = x1(i) - x2(j);
        x_right = dx + a;  
        x_left = dx - a;
        Z1_to_Z2(i,j) = coef * K2_function(a, x_right, h, x_left, h);
    end
end

% symmetry
Z2_to_Z1 = Z1_to_Z2.';          

% Make Z matrix
Z = [Z1_self, Z1_to_Z2; Z2_to_Z1, Z2_self];
Z = Z - Z(1,:);
Z(1,:) = delta;

% Make V matrix
V = delta * [zeros(N1,1) ; (v1-v2)*ones(N2,1) ];

% Surface charge density
rho = Z \ V;
rho1 = rho(1:N1); 
rho2 = rho(N1+1:end);

Q_prime = sum(rho1) * delta;
C_prime = Q_prime / (v1 - v2);

% plotting rho
figure;
plot(x1, rho1, 'LineWidth', 1.2); hold on;
plot(x2, rho2, 'LineWidth', 1.2); grid on;
xlabel('x position (m)'); ylabel('\rho_s (C/m^2)');
legend('Conductor (y = h)', 'Ground (y = 0)', 'Location','best');
title(sprintf('Galerkin (Pulse–Pulse) — w_1=%.3g, w_2=%.3g, h=%.3g, N_1=%d, N_2=%d', ...
      w1, w2, h, N1, N2));

% several characteristic values of w1, w2, and h
w1_list = [0.5, 1.0];
w2_list = [2.0, 10.0];
h_list  = [0.5, 1.0, 2.0];

% --- MATCHED TO POINT-MATCHING N1_list ---
N_list  = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300];

% results: [w1, w2_used, h, N1, N2, Ntot, Cprime_MoM, Cprime_emp, Err%]
max_rows = length(w1_list)*length(w2_list)*length(h_list)*length(N_list);
results = zeros(max_rows,9);
row = 1;

for i_a = 1:length(w1_list)
  for i_b = 1:length(w2_list)
    for i_c = 1:length(h_list)

      w1_cur = w1_list(i_a);
      w2_cur = w2_list(i_b);
      h_cur  = h_list(i_c);

      ratio = w1_cur / h_cur;
      if ratio < 1
        p = 0.04*(1 - ratio)^2;
      else
        p = 0;
      end
      eps_eff = (epsilon_r + 1)/2 + ((epsilon_r - 1)/2)*((1 + 12/ratio)^(-1/2) + p);
      if ratio <= 1
        Z0_emp = eta0/(2*pi*sqrt(eps_eff)) * log(8/ratio + ratio/4);
      else
        Z0_emp = eta0/sqrt(eps_eff) * (ratio + 1.393 + 0.667 * log(ratio + 1.444))^-1;
      end
      Cprime_emp = 1/(Z0_emp*c0);

      for i_n = 1:length(N_list)

        N1_cur = N_list(i_n);
        delta = w1_cur / N1_cur;
        a = delta/2;
        N2_cur = max(1, round(w2_cur / delta));
        w2_used = N2_cur * delta;

        % Re-establish X1 and and X2 coordinates with new variables
        x1 = linspace(-w1_cur/2 + a, +w1_cur/2 - a, N1_cur).';
        x2 = linspace(-w2_used/2 + a, +w2_used/2 - a, N2_cur).';

        Z1_self = zeros(N1_cur,N1_cur);
        for i = 1:N1_cur
          for j = 1:N1_cur
            dx = x1(i) - x1(j);
            x_right = dx + a;  
            x_left = dx - a;
            Z1_self(i,j) = coef * K2_function(a, x_right, 0, x_left, 0);
          end
        end

        Z2_self = zeros(N2_cur,N2_cur);
        for i = 1:N2_cur
          for j = 1:N2_cur
            dx = x2(i) - x2(j);
            x_right = dx + a;  
            x_left = dx - a;
            Z2_self(i,j) = coef * K2_function(a, x_right, 0, x_left, 0);
          end
        end

        Z1_to_Z2 = zeros(N1_cur,N2_cur);
        for i = 1:N1_cur
          for j = 1:N2_cur
            dx = x1(i) - x2(j);
            x_right = dx + a;  
            x_left = dx - a;
            Z1_to_Z2(i,j) = coef * K2_function(a, x_right, h_cur, x_left, h_cur);
          end
        end

        % symmetry
        Z2_to_Z1 = Z1_to_Z2.';  

        % Make z matrix
        Z = [Z1_self, Z1_to_Z2; Z2_to_Z1, Z2_self];
        Z = Z - Z(1,:);
        Z(1,:) = delta;

        % make V matrix
        V = delta * [zeros(N1_cur,1); (v1-v2) * ones(N2_cur,1)];

        % solve, top charge -> C'
        rho = Z \ V;
        rho1 = rho(1:N1_cur);
        Q_prime = sum(rho1) * delta;
        C_prime_Galerkin = Q_prime/(v1 - v2);

        % percent error
        percent_error = abs(C_prime_Galerkin - Cprime_emp)/Cprime_emp*100;

        % store row
        results(row,1) = w1_cur;
        results(row,2) = w2_used;
        results(row,3) = h_cur;
        results(row,4) = N1_cur;
        results(row,5) = N2_cur;
        results(row,6) = N1_cur + N2_cur;
        results(row,7) = C_prime_Galerkin;
        results(row,8) = Cprime_emp;
        results(row,9) = percent_error;
        row = row + 1;
      end
    end
  end
end

results = results(1:row-1,:);

fprintf('\nCapacitance and Relative Error (MoM vs Empirical)\n');
fprintf(' w1     w2_used   h      N1    N2   Ntot        C''_MoM(F/m)      C''_emp(F/m)    Err(%%)\n');
for r = 1:size(results,1)
  fprintf('%5.2f   %7.2f  %5.2f  %5d  %5d  %5d   %14.6e  %14.6e   %7.3f\n', ...
    results(r,1), results(r,2), results(r,3), ...
    results(r,4), results(r,5), results(r,6), ...
    results(r,7), results(r,8), results(r,9));
end

% Plot: empirical % error vs N1
figure; hold on; grid on;
labels = {};
idx = 1;
for i_a = 1:length(w1_list)
  for i_b = 1:length(w2_list)
    for i_c = 1:length(h_list)
      Ns = zeros(length(N_list),1);
      Es = zeros(length(N_list),1);
      for k = 1:length(N_list)
        Nwant = N_list(k);
        for r = 1:size(results,1)
          if results(r,1)==w1_list(i_a) && results(r,3)==h_list(i_c) && results(r,4)==Nwant
            Ns(k) = results(r,4);
            Es(k) = results(r,9);
          end
        end
      end
      plot(Ns, Es, '-o', 'LineWidth', 1.2);
      labels{idx} = sprintf('w1=%.2f, w2=%.2f, h=%.2f', w1_list(i_a), w2_list(i_b), h_list(i_c));
      idx = idx + 1;
    end
  end
end
xlabel('N_1'); ylabel('Relative % Error vs empirical');
title('Relative % Error (empirical) vs N_1');
legend(labels, 'Location', 'best'); hold off;

% Convergence graph: relative error to C_inf (log–log vs N1)
triplets = unique(results(:,[1,2,3]), 'rows', 'stable');  % [w1, w2_used, h]
figure; hold on; grid on;
leg = strings(size(triplets,1),1);
for t = 1:size(triplets,1)
    w1c = triplets(t,1);  w2c = triplets(t,2);  hc = triplets(t,3);
    mask = (results(:,1)==w1c) & (results(:,2)==w2c) & (results(:,3)==hc);
    sub  = results(mask, :);
    [Nvec, idx] = sort(sub(:,4));
    Cvec = sub(idx,7);   % C'_MoM
    if numel(Nvec) < 2, continue; end
    N1a = Nvec(end-1);  N1b = Nvec(end);
    C1a = Cvec(end-1);  C1b = Cvec(end);
    Cinf = (N1b*C1b - N1a*C1a) / (N1b - N1a);
    relErr = abs(Cvec - Cinf) ./ abs(Cinf);
    p = polyfit(log(Nvec), log(relErr), 1);
    slope_est = p(1);
    loglog(Nvec, relErr, '-o', 'LineWidth', 1.2);
    leg(t) = sprintf('w1=%.2g, w2=%.2g, h=%.2g | C_\\infty=%.3e, slope≈%.2f', ...
                     w1c, w2c, hc, Cinf, slope_est);
end
xlabel('N_1');
ylabel('Relative error to C''_\infty');
title('Convergence of MoM Capacitance (relative error to C_\infty, log–log)');
legend(leg(leg~=""), 'Location','southwest');
hold off;

% 1/N reference
ax = gca;
xGuide = get(ax, 'XLim');
xg = logspace(log10(xGuide(1)), log10(xGuide(2)), 50);
yGuide = get(ax, 'YLim');
k = 0.5 * yGuide(2);
hold on;
loglog(xg, k./xg, '--', 'LineWidth', 1.0);
text(exp(mean(log(xGuide))), k/exp(mean(log(xGuide))), '  ~1/N', 'VerticalAlignment','bottom');
hold off;

%% ---------- Galerkin: convergence RATE vs Delta (assignment formula) ----------
% r_i = log2(err_i/err_{i+1}) / log2(Δ_i/Δ_{i+1}),  err_i = |C'(N1_i) - C_inf|,  Δ = w1/N1
% Group by (w1,h) and sub-group by rounded w2_used to avoid mixing grounds.

triplets_wh = unique(results(:,[1,3]), 'rows', 'stable');  % [w1, h]

figure('Name','Galerkin: convergence rate vs Delta'); hold on; grid on;
xlabel('\Delta = w_1/N_1 (m)'); set(gca,'XScale','log'); set(gca,'XDir','reverse'); % coarse→fine left→right
ylabel('Convergence rate  r_i');
title('Galerkin convergence rate  r_i = log_2(err_i/err_{i+1}) / log_2(\Delta_i/\Delta_{i+1})');

for t = 1:size(triplets_wh,1)
    w1c = triplets_wh(t,1);  hc = triplets_wh(t,2);

    mask_wh = (results(:,1)==w1c) & (results(:,3)==hc);
    sub_wh  = results(mask_wh,:);  % [w1 w2_used h N1 N2 Ntot C_MoM C_emp Err%]
    if isempty(sub_wh), continue; end

    % sub-group by rounded w2_used (meters)
    w2_bins = unique(round(sub_wh(:,2), 6));
    for b = 1:numel(w2_bins)
        w2b = w2_bins(b);
        sub  = sub_wh( round(sub_wh(:,2),6) == w2b, : );
        if size(sub,1) < 3, continue; end

        [Nvec, idx] = sort(sub(:,4));   % N1 ascending
        Cvec = sub(idx,7);
        Dlt  = w1c ./ Nvec;             % Δ = w1/N1

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
            fprintf('\n[Galerkin] (w1=%.2g, h=%.2g; w2≈%.6g)\n', w1c, hc, w2b);
            fprintf('  i    N1_i   N1_{i+1}        Δ_i          Δ_{i+1}      r_i\n');
            for i=1:numPairs
                fprintf('%3d  %6d   %8d   %12.5e  %12.5e   %7.3f\n', ...
                    i, Nvec(idd(i)), Nvec(idd(i+1)), Ddesc(i), Ddesc(i+1), rates(i));
            end
        end
    end
end
legend('show','Location','bestoutside'); hold off;
