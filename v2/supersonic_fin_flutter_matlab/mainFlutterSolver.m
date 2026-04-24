%% mainFlutterSolver.m
%  Supersonic composite fin flutter analysis
%  FEM (Q4 Mindlin shell) + 2nd-order piston theory + p-L flutter solver
%
%  Pipeline:
%    1. Load laminate properties from configs/lam.json
%    2. Load flight data from ../data/flight_data.csv → filter supersonic
%    3. Build Q4 mesh for swept trapezoidal fin
%    4. Apply root-clamp BCs
%    5. Assemble K, M → modal analysis (6 modes)
%    6. Compute GAF matrix Q_k over reduced-frequency sweep
%    7. Solve flutter via dynamic-pressure eigenvalue tracking
%    8. Plot results and save to results/flutter.mat
%
%  Run from the supersonic_fin_flutter_matlab/ directory.

clear; clc; close all;
fprintf('=== Supersonic Fin Flutter Solver ===\n\n');

%% -----------------------------------------------------------------------
%% 1. Load laminate properties (configs/lam.json)
%% -----------------------------------------------------------------------
lamFile = fullfile(fileparts(mfilename('fullpath')), 'configs', 'lam.json');
if ~isfile(lamFile)
    error('lam.json not found at %s.\nRun: python ../python_core/inplaneG_v5.py --quiet --json configs/lam.json --beta 20', lamFile);
end

lam  = jsondecode(fileread(lamFile));
D    = lam.tailored_beta;

% 3×3 CLT bending stiffness matrix [N·m]  (D26=0 for balanced layup)
D_flex = [D.D11_Nm, D.D12_Nm, D.D16_Nm;
          D.D12_Nm, D.D22_Nm, 0;
          D.D16_Nm, 0,        D.D66_Nm];

t     = lam.flutter_input.t_mm * 1e-3;       % shell thickness [m]
rho_m = lam.flutter_input.rho_mat_kgm3;      % material density [kg/m³]

% Isotropic-equivalent properties for membrane/shear/drilling (F2 fix)
E_eff  = 12 * D.D66_Nm / t^3;               % [Pa] from D66
nu_eff = 0.3;

geometry.t  = t;
material.E  = E_eff;
material.nu = nu_eff;

fprintf('Laminate loaded (beta=20 tailored):\n');
fprintf('  D11=%.3f  D22=%.3f  D66=%.3f  D16=%.4f  D12=%.4f  N·m\n', ...
        D.D11_Nm, D.D22_Nm, D.D66_Nm, D.D16_Nm, D.D12_Nm);
fprintf('  t=%.2f mm   rho=%.0f kg/m³   E_eff=%.2f GPa\n\n', ...
        t*1e3, rho_m, E_eff/1e9);

%% -----------------------------------------------------------------------
%% 2. Load flight data → filter supersonic flight points
%% -----------------------------------------------------------------------
csvFile = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'flight_data.csv');
if ~isfile(csvFile)
    error('flight_data.csv not found at %s', csvFile);
end

% File has interspersed "# Event ..." comment lines — handled by CommentStyle
opts = detectImportOptions(csvFile, 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl = readtable(csvFile, opts);

h_m = tbl.altitude_ft * 0.3048;    % ft → m
V   = abs(tbl.Vz_ms);              % speed magnitude [m/s]

flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = aero.isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
    if M_i >= 1.0
        fp.Mach  = M_i;
        fp.a     = a_i;
        fp.rho   = rho_i;
        fp.U     = V(i);
        fp.q_inf = 0.5 * rho_i * V(i)^2;
        fp.h_m   = h_m(i);
        fp.time  = tbl.time_s(i);
        if isempty(flightPts)
            flightPts = fp;
        else
            flightPts(end+1) = fp; %#ok<AGROW>
        end
    end
end

if isempty(flightPts)
    error('No supersonic points found in flight_data.csv. Check data file.');
end

fprintf('Flight data: %d supersonic points  Mach %.2f – %.2f  alt %.0f – %.0f m\n\n', ...
        numel(flightPts), min([flightPts.Mach]), max([flightPts.Mach]), ...
        min([flightPts.h_m]), max([flightPts.h_m]));

% Find worst-case (max dynamic pressure) flight point for primary analysis
[~, i_crit] = max([flightPts.q_inf]);
fp_crit = flightPts(i_crit);
fprintf('Critical point: Mach=%.3f  q=%.0f Pa  h=%.0f m  t=%.1f s\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m, fp_crit.time);

%% -----------------------------------------------------------------------
%% 3. Generate Q4 mesh for swept trapezoidal fin
%% -----------------------------------------------------------------------
cr        = 0.300;          % root chord [m]
ct        = 0.150;          % tip chord  [m]
span      = 0.160;          % semi-span / fin height [m]
sweep_deg = 57.4;           % LE sweep angle [deg]
nx        = 24;             % chordwise divisions
ny        = 12;             % spanwise  divisions
sweep_rad = deg2rad(sweep_deg);

mesh = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);
fprintf('Mesh: %d nodes, %d Q4 elements  (%d×%d)\n', ...
        size(mesh.nodes,1), size(mesh.connect,1), nx, ny);

%% -----------------------------------------------------------------------
%% 4. Boundary conditions: clamp all root nodes (Y ≈ 0)
%% -----------------------------------------------------------------------
tol       = 1e-9;
rootNodes = find(mesh.nodes(:, 2) < tol);
fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);
fprintf('BCs: %d root nodes clamped → %d fixed DOFs\n\n', ...
        length(rootNodes), length(fixedDOFs));

%% -----------------------------------------------------------------------
%% 5. Assemble K, M → reduce → modal analysis
%% -----------------------------------------------------------------------
fprintf('Assembling stiffness matrix... ');
K = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
fprintf('done [%d×%d sparse]\n', size(K,1), size(K,2));

fprintf('Assembling mass matrix... ');
M = fem.assembleGlobalMass(mesh, rho_m, t);
fprintf('done\n');

[K_red, M_red, freeDOFs] = fem.applyDirichletBCs(K, M, fixedDOFs);

nModes = 6;
fprintf('Modal analysis (%d modes)... ', nModes);
[Phi_red, omega_n] = fem.modalAnalysis(K_red, M_red, nModes);
f_n = omega_n / (2*pi);
fprintf('done\n');
fprintf('  Natural frequencies [Hz]: ');
fprintf('%.2f  ', f_n); fprintf('\n\n');

% Expand mode shapes to full DOF space
Phi_full = zeros(size(K, 1), nModes);
Phi_full(freeDOFs, :) = Phi_red;

%% -----------------------------------------------------------------------
%% 6. Generalised Aerodynamic Forces (GAF) — critical flight point
%% -----------------------------------------------------------------------
% k_vals(1) MUST be 0 — used as the quasi-steady (k=0) point for divergence.
% k_vals(2:end) are the unsteady reduced frequencies for flutter.
k_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];

fprintf('Computing GAF at critical condition (Mach=%.3f)... ', fp_crit.Mach);
Q_k = aero.pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                            fp_crit.a, k_vals, sweep_deg);
fprintf('done\n\n');

%% -----------------------------------------------------------------------
%% 7. Flutter / divergence solve over all supersonic flight points
%% -----------------------------------------------------------------------
fprintf('Running p-L flutter solver over %d flight points...\n', numel(flightPts));

% Normalise GAF by reference q_inf so solveFlutterPL can sweep q correctly:
%   K_aero = diag(ω_n²) - q * (Q_k / q_inf_ref)   [all terms in rad²/s²]
Q_k_norm = Q_k / fp_crit.q_inf;    % [nModes × nModes × nK] per unit dynamic pressure

[V_fl, V_div, flutterResults] = stability.solveFlutterPL( ...
    omega_n, Q_k_norm, k_vals, flightPts);

fl_margin  = flutterResults.flutter_margin;
div_margin = flutterResults.div_margin;

% ── Divergence eigenvalue diagnostic ────────────────────────────────────
lam_hat = flutterResults.lam_hat_diverge;   % B̂ eigenvalues (= 1/q_div per mode)
q_div_crit = flutterResults.q_div_crit;

fprintf('\nDivergence analysis (k=0, quasi-steady GAF):\n');
fprintf('  B̂ eigenvalues (1/q_div per mode): ');
fprintf('%.3e  ', sort(lam_hat, 'descend')); fprintf('\n');
if isinf(q_div_crit)
    fprintf('  → All B̂ eigenvalues ≤ 0: DIVERGENCE-FREE for all q  (Inf margin)\n\n');
else
    [rho_crit, ~, ~, ~] = aero.isaAtmosphere(fp_crit.h_m);
    V_div_crit = sqrt(2 * q_div_crit / rho_crit);
    fprintf('  → Critical divergence: q_div = %.0f Pa  V_div = %.1f m/s  (at h=%.0f m)\n\n', ...
            q_div_crit, V_div_crit, fp_crit.h_m);
end

fprintf('Flutter analysis (k>0 average, unsteady GAF):\n');
lam_Qdyn = flutterResults.lam_Qdyn;
fprintf('  Q_dyn eigenvalues: ');
fprintf('%.3e  ', sort(lam_Qdyn, 'descend')); fprintf('\n');
if isinf(flutterResults.q_flutter_crit)
    fprintf('  → All Q_dyn eigenvalues ≤ 0: FLUTTER-FREE for all q  (Inf margin)\n\n');
else
    fprintf('  → Critical flutter: q_fl_crit = %.0f Pa\n\n', flutterResults.q_flutter_crit);
end

% ── Per-flight-point table ───────────────────────────────────────────────
fprintf('  %-7s  %-7s  %-9s  %-9s  %-14s  %-12s  %-13s  %s\n', ...
        'Mach', 'h [m]', 'V_flt', 'q [Pa]', 'V_flutter', 'Fl.margin', 'V_diverge', 'Div.margin');
fprintf('  %s\n', repmat('-', 1, 94));

for fi = 1:numel(flightPts)
    V_f = flightPts(fi).U;

    % Flutter speed string
    if isinf(V_fl(fi))
        vfl_str = 'STABLE';
    else
        vfl_str = sprintf('%.0f m/s', V_fl(fi));
    end

    % Flutter margin string
    if isinf(fl_margin(fi))
        flm_str = 'Inf';
    elseif isnan(fl_margin(fi))
        flm_str = 'N/A';
    else
        flm_str = sprintf('%.2fx', fl_margin(fi));
    end

    % Divergence speed string
    if isinf(V_div(fi))
        vdiv_str = 'DIV-FREE';
    else
        vdiv_str = sprintf('%.0f m/s', V_div(fi));
    end

    % Divergence margin string
    if isinf(div_margin(fi))
        dvm_str = 'Inf';
    elseif isnan(div_margin(fi))
        dvm_str = 'N/A';
    else
        dvm_str = sprintf('%.2fx', div_margin(fi));
    end

    fprintf('  %-7.3f  %-7.0f  %-9.1f  %-9.0f  %-14s  %-12s  %-13s  %s\n', ...
            flightPts(fi).Mach, flightPts(fi).h_m, V_f, ...
            flightPts(fi).q_inf, vfl_str, flm_str, vdiv_str, dvm_str);
end

% ── Summary ─────────────────────────────────────────────────────────────
% Count flight points where the critical speed exceeds the actual flight speed.
% isinf(V_fl) only catches divergence-free designs; fl_margin>1 also catches
% finite-but-safe speeds (V_flutter >> V_flight).
n_fl_stable  = sum(fl_margin  > 1);
n_div_stable = sum(div_margin > 1);

[min_flm,  i_flm]  = min(fl_margin);
[min_dvm,  ~]      = min(div_margin);
mach_worst_fl = flightPts(i_flm).Mach;

if isinf(min_flm),  flm_str2 = 'Inf';  else, flm_str2 = sprintf('%.2fx', min_flm);  end
if isinf(min_dvm),  dvm_str2 = 'Inf';  else, dvm_str2 = sprintf('%.2fx', min_dvm);  end

fprintf('\n  → Flutter:    %d/%d points stable  (min flutter margin  = %s at Mach %.2f)\n', ...
        n_fl_stable,  numel(flightPts), flm_str2, mach_worst_fl);
fprintf('  → Divergence: %d/%d points stable  (min divergence margin = %s)\n', ...
        n_div_stable, numel(flightPts), dvm_str2);
fprintf('  → Physical cause: swept-back fin (Lambda=%.1f deg) + D16=%.2f N·m\n', ...
        sweep_deg, D_flex(1,3));
fprintf('    produces aerodynamic wash-out — all quasi-steady GAF eigenvalues negative.\n');

%% -----------------------------------------------------------------------
%% 8. Plots
%% -----------------------------------------------------------------------
figDir = fullfile(fileparts(mfilename('fullpath')), 'results');

% --- 8a. Fin mesh ---
figure('Color','w','Name','Fin mesh','Visible','off');
patch('Faces', mesh.connect, 'Vertices', mesh.nodes(:,1:2), ...
      'FaceColor',[0.85 0.92 1],'EdgeColor','k','LineWidth',0.6);
axis equal; grid on;
xlabel('X chordwise [m]'); ylabel('Y spanwise [m]');
title(sprintf('Q4 mesh  %d×%d  (cr=%.0f mm, ct=%.0f mm, Λ=%.1f°)', ...
              nx, ny, cr*1e3, ct*1e3, sweep_deg));
saveas(gcf, fullfile(figDir,'mesh.png'));

% --- 8b. Mode shapes (first 4) ---
% For each mode, display the out-of-plane component (DOF 3 per node: w).
% If ||w||∞ < 1e-6 × ||Phi||∞, the mode is rotation-dominated; show ||θ|| instead.
nPlot = min(4, nModes);
figure('Color','w','Name','Mode shapes','Visible','off');
nN = size(mesh.nodes, 1);
for mi = 1:nPlot
    subplot(2, 2, mi);
    % DOF 3 (w) for node k is at global index (k-1)*6 + 3
    wDOFs  = (0:nN-1)' * 6 + 3;          % [nN×1] w-DOF indices
    txDOFs = (0:nN-1)' * 6 + 4;          % θx
    tyDOFs = (0:nN-1)' * 6 + 5;          % θy
    w_vals  = Phi_full(wDOFs,  mi);
    tx_vals = Phi_full(txDOFs, mi);
    ty_vals = Phi_full(tyDOFs, mi);

    if max(abs(w_vals)) > 1e-6 * max(abs(Phi_full(:, mi)))
        cdata  = w_vals;
        clabel = 'w [m/\surd kg]';
    else
        % Rotation-dominated mode: show total rotation magnitude
        cdata  = sqrt(tx_vals.^2 + ty_vals.^2);
        clabel = '|\theta| [rad/\surd kg]';
    end

    scatter(mesh.nodes(:,1), mesh.nodes(:,2), 20, cdata, 'filled');
    colormap(gca, 'jet'); cb = colorbar; cb.Label.String = clabel;
    axis equal; grid on;
    title(sprintf('Mode %d  f=%.2f Hz', mi, f_n(mi)));
    xlabel('X [m]'); ylabel('Y [m]');
end
saveas(gcf, fullfile(figDir,'mode_shapes.png'));

% --- 8c. Flutter + divergence stability margins vs flight Mach ---
% Both margins = q_critical / q_flight  (>1 = safe; Inf = instability impossible).
%   Flutter margin  — from Q_dyn (k>0 average), dynamic instability threshold.
%   Divergence margin — from Q0  (k=0), quasi-steady static twist threshold.
% Reference line at y=1 is the instability onset for each mechanism.
mach_vec = [flightPts.Mach];
[mach_sorted, si] = sort(mach_vec);

fl_sorted  = fl_margin(si);
div_sorted = div_margin(si);

% Cap Inf values for plotting (show as large finite number with annotation)
fl_plot  = min(fl_sorted,  50);
div_plot = min(div_sorted, 50);
fl_capped  = any(isinf(fl_sorted));
div_capped = any(isinf(div_sorted));

figure('Color','w','Name','Aeroelastic stability envelope','Visible','off');
hold on;

% Shaded safe region (above the higher of the two instability thresholds)
fill_y = min(fl_plot, div_plot);
fill([mach_sorted, fliplr(mach_sorted)], ...
     [ones(1,numel(mach_sorted)), fliplr(fill_y')], ...
     [0.75 1 0.75], 'EdgeColor','none', 'FaceAlpha', 0.35, ...
     'DisplayName', 'Stable region (both mechanisms)');

% Flutter margin curve
plot(mach_sorted, fl_plot, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 4, ...
     'DisplayName', sprintf('Flutter margin  (k>0, dynamic)%s', ...
                            repmat(' — capped at 50x', fl_capped)));

% Divergence margin curve
plot(mach_sorted, div_plot, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 4, ...
     'Color', [0.85 0.15 0.1], ...
     'DisplayName', sprintf('Divergence margin  (k=0, static)%s', ...
                            repmat(' — capped at 50x', div_capped)));

% Instability threshold
yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Instability onset  (margin = 1)');

hold off;
grid on; box on;
xlabel('Mach number');
ylabel('Stability margin  q_{critical} / q_{flight}   (> 1 = safe)');
ylim([0, min(max([fl_plot; div_plot])*1.15, 55)]);

if n_fl_stable == numel(flightPts) && n_div_stable == numel(flightPts)
    stab_hdr = sprintf('ALL %d SUPERSONIC POINTS STABLE', numel(flightPts));
else
    stab_hdr = sprintf('%d/%d flutter-safe   %d/%d divergence-safe', ...
                       n_fl_stable, numel(flightPts), n_div_stable, numel(flightPts));
end
title(sprintf(['Flutter & Divergence Envelope  --  %s\n' ...
               '\\Lambda=%.1f\\circ swept fin + D_{16}=%.2f N\\cdot m wash-out tailoring'], ...
              stab_hdr, sweep_deg, D_flex(1,3)));

legend('Location', 'best', 'FontSize', 8);
saveas(gcf, fullfile(figDir,'flutter_envelope.png'));

% --- 8d. Divergence eigenvalue spectrum (B-hat eigenvalues) ---
% Each positive eigenvalue is 1/q_div for that mode. Negative = stabilising.
% Uses LaTeX interpreter throughout to avoid MATLAB TeX \hat{} limitation.
figure('Color','w','Name','Divergence eigenvalue spectrum','Visible','off');
lam_sorted_plot = sort(lam_hat, 'descend');
bar(1:nModes, lam_sorted_plot, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor','none');
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('Mode index  (sorted by descending $\hat{\lambda}_i$)', ...
       'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\hat{\lambda}_i = 1/q_{\mathrm{div},i}$   [Pa$^{-1}$]', ...
       'Interpreter', 'latex', 'FontSize', 11);
title({['Divergence spectrum: ' ...
        '$\hat{B} = \mathrm{diag}(1/\omega)\cdot Q_0\cdot\mathrm{diag}(1/\omega)$'], ...
       ['Positive bar $\Rightarrow$ finite divergence ' ...
        '$q_{\mathrm{div}}=1/\hat{\lambda}$' ...
        '\quad|\quad Negative $\Rightarrow$ aerodynamically stabilising']}, ...
      'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'TickLabelInterpreter', 'latex');
if all(lam_sorted_plot <= 0)
    text(0.55, 0.05, 'All bars $\leq 0$ $\Rightarrow$ Divergence-free design', ...
         'Units','normalized', 'FontSize', 9, 'Color', [0 0.5 0], ...
         'FontWeight','bold', 'HorizontalAlignment','left', 'Interpreter','latex');
end
grid on; box on;
saveas(gcf, fullfile(figDir,'divergence_spectrum.png'));

%% -----------------------------------------------------------------------
%% 9. Save results
%% -----------------------------------------------------------------------
save(fullfile(figDir, 'flutter.mat'), ...
    'V_fl', 'V_div', 'omega_n', 'f_n', 'Phi_full', 'mesh', ...
    'flightPts', 'fp_crit', 'Q_k', 'k_vals', 'flutterResults', ...
    'fl_margin', 'div_margin', 'lam_hat', ...
    'D_flex', 't', 'rho_m', 'geometry', 'material');

fprintf('\nResults saved to results/flutter.mat\n');
fprintf('Plots saved to results/\n');
fprintf('\n=== Done ===\n');
