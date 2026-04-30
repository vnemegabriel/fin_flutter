%% mainFlutterSolver.m
%  Supersonic composite fin flutter analysis
%  FEM (Q4 Mindlin shell) + 2nd-order piston theory + p-L flutter solver
%
%  Pipeline:
%    1. Build laminate properties via CLT + Halpin-Tsai (+core/buildCLTLaminate.m)
%    2. Load flight data from data/flight_data.csv → filter supersonic points
%    3. Build Q4 mesh for swept trapezoidal fin
%    4. Apply root-clamp BCs
%    5. Assemble K, M → modal analysis (6 modes)
%    6. Compute GAF matrix Q_k over reduced-frequency sweep
%    7. Solve flutter/divergence via dynamic-pressure eigenvalue tracking
%    8. Plot results and save to results/flutter.mat
%
%  Run from the supersonic_fin_flutter_matlab/ directory.

clear; clc; close all;
fprintf('=== Supersonic Fin Flutter Solver ===\n\n');

%% -----------------------------------------------------------------------
%% 1. Build laminate (AR1 T700/Epoxy, beta=20° aeroelastic tailoring)
%% -----------------------------------------------------------------------
% CLT bending stiffness + Halpin-Tsai micromechanics, no external file needed.
% beta=20° rotates all plies globally, introducing D16 ≠ 0 (bend-twist coupling).
% The resulting wash-out mechanism (bending reduces local AoA) stabilises flutter.
[D_flex, t, rho_lam, D_info] = core.buildCLTLaminate(20, 0.50);

% Isotropic-equivalent stiffness for Q4 membrane/shear/drilling DOFs.
% Using D66 gives an effective shear modulus consistent with CLT torsional stiffness.
E_eff  = 12 * D_flex(3,3) / t^3;   % [Pa]
nu_eff = 0.3;

geometry.t  = t;
material.E  = E_eff;
material.nu = nu_eff;

fprintf('Laminate (AR1 T700/Epoxy, beta=20 tailored):\n');
fprintf('  D11=%.3f  D22=%.3f  D66=%.3f  D16=%.4f  D12=%.4f  N·m\n', ...
        D_info.D11_Nm, D_info.D22_Nm, D_info.D66_Nm, D_info.D16_Nm, D_info.D12_Nm);
fprintf('  t=%.2f mm   rho=%.0f kg/m³   E_eff=%.2f GPa\n\n', ...
        t*1e3, rho_lam, E_eff/1e9);

%% -----------------------------------------------------------------------
%% 2. Load flight data → filter supersonic flight points
%% -----------------------------------------------------------------------
csvFile = fullfile(fileparts(mfilename('fullpath')), 'data', 'flight_data.csv');
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
    % M >= 1.05 ensures beta = sqrt(M^2-1) >= 0.312; piston theory validity per Lighthill (1953) and NACA TN 4021
    if M_i >= 1.05
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

fprintf('Flight data: %d supersonic points (M >= 1.05)  Mach %.2f – %.2f  alt %.0f – %.0f m\n\n', ...
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
M = fem.assembleGlobalMass(mesh, rho_lam, t);
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

% Normalise GAF by reference q_inf so solveFlutterPL can sweep q correctly
Q_k_norm = Q_k / fp_crit.q_inf;    % [nModes × nModes × nK] per unit dynamic pressure

[V_fl, V_div, flutterResults] = stability.solveFlutterPL( ...
    omega_n, Q_k_norm, k_vals, flightPts);

fl_margin  = flutterResults.flutter_margin;
div_margin = flutterResults.div_margin;

% ── Stability summary at critical flight point ───────────────────────────
q_div_crit = flutterResults.q_div_crit;

fprintf('\nStability at critical flight point (Mach=%.3f  q=%.0f Pa):\n', ...
        fp_crit.Mach, fp_crit.q_inf);

% Divergence (quasi-steady, k=0):
%   Checks whether aerodynamic twist moments ever overcome structural stiffness
%   for a static load increase. A negative divergence eigenvalue means the
%   aerodynamic load acts as a restoring moment — the fin twists to reduce lift.
if isinf(q_div_crit)
    fprintf('  Divergence:  STABLE at all speeds\n');
    fprintf('    Swept-back geometry (%.1f deg) + D16=%.2f N·m bend-twist coupling\n', ...
            sweep_deg, D_flex(1,3));
    fprintf('    produce aero wash-out: bending always decreases local AoA.\n\n');
else
    [rho_crit, ~, ~, ~] = aero.isaAtmosphere(fp_crit.h_m);
    V_div_crit = sqrt(2 * q_div_crit / rho_crit);
    fprintf('  Divergence:  q_div = %.0f Pa   V_div = %.1f m/s  (at h=%.0f m)\n\n', ...
            q_div_crit, V_div_crit, fp_crit.h_m);
end

% Flutter (unsteady, k>0):
%   Checks whether coupled aero-structural oscillations grow unboundedly.
%   Stability requires the unsteady coupling matrix to have no energy-adding eigenvalues.
if isinf(flutterResults.q_flutter_crit)
    fprintf('  Flutter:     STABLE at all speeds\n');
    fprintf('    All unsteady aerodynamic coupling terms are stabilising.\n\n');
else
    fprintf('  Flutter:     q_fl_crit = %.0f Pa\n\n', flutterResults.q_flutter_crit);
end

% ── Per-flight-point table ───────────────────────────────────────────────
fprintf('  %-7s  %-7s  %-9s  %-9s  %-14s  %-12s  %-13s  %s\n', ...
        'Mach', 'h [m]', 'V_flt', 'q [Pa]', 'V_flutter', 'Fl.margin', 'V_diverge', 'Div.margin');
fprintf('  %s\n', repmat('-', 1, 94));

for fi = 1:numel(flightPts)
    V_f = flightPts(fi).U;

    if isinf(V_fl(fi))
        vfl_str = 'STABLE';
    else
        vfl_str = sprintf('%.0f m/s', V_fl(fi));
    end

    if isinf(fl_margin(fi))
        flm_str = 'Inf';
    elseif isnan(fl_margin(fi))
        flm_str = 'N/A';
    else
        flm_str = sprintf('%.2fx', fl_margin(fi));
    end

    if isinf(V_div(fi))
        vdiv_str = 'DIV-FREE';
    else
        vdiv_str = sprintf('%.0f m/s', V_div(fi));
    end

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
n_fl_stable  = sum(fl_margin  > 1);
n_div_stable = sum(div_margin > 1);

[min_flm,  i_flm]  = min(fl_margin);
[min_dvm,  ~]      = min(div_margin);
mach_worst_fl = flightPts(i_flm).Mach;

if isinf(min_flm),  flm_str2 = 'Inf';  else, flm_str2 = sprintf('%.2fx', min_flm);  end
if isinf(min_dvm),  dvm_str2 = 'Inf';  else, dvm_str2 = sprintf('%.2fx', min_dvm);  end

fprintf('\n  → Flutter:    %d/%d points stable  (min margin = %s at Mach %.2f)\n', ...
        n_fl_stable,  numel(flightPts), flm_str2, mach_worst_fl);
fprintf('  → Divergence: %d/%d points stable  (min margin = %s)\n\n', ...
        n_div_stable, numel(flightPts), dvm_str2);

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
% Displays the out-of-plane displacement (w) for bending-dominated modes,
% or rotation magnitude for torsion-dominated modes.
nPlot = min(4, nModes);
figure('Color','w','Name','Mode shapes','Visible','off');
nN = size(mesh.nodes, 1);
for mi = 1:nPlot
    subplot(2, 2, mi);
    wDOFs  = (0:nN-1)' * 6 + 3;
    txDOFs = (0:nN-1)' * 6 + 4;
    tyDOFs = (0:nN-1)' * 6 + 5;
    w_vals  = Phi_full(wDOFs,  mi);
    tx_vals = Phi_full(txDOFs, mi);
    ty_vals = Phi_full(tyDOFs, mi);

    if max(abs(w_vals)) > 1e-6 * max(abs(Phi_full(:, mi)))
        cdata  = w_vals;
        clabel = 'w [m/\surd kg]';
    else
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
% Stability margin = q_critical / q_flight (> 1 = safe, Inf = unconditionally stable).
% A margin of 2 means the fin could withstand twice the actual flight dynamic pressure
% before going unstable. Both flutter (dynamic oscillation) and divergence (static
% twist-up) are checked; the more restrictive determines the design limit.
mach_vec = [flightPts.Mach];
[mach_sorted, si] = sort(mach_vec);

fl_sorted  = fl_margin(si);
div_sorted = div_margin(si);

fl_plot  = min(fl_sorted,  50);
div_plot = min(div_sorted, 50);
fl_capped  = any(isinf(fl_sorted));
div_capped = any(isinf(div_sorted));

figure('Color','w','Name','Aeroelastic stability envelope','Visible','off');
hold on;

fill_y = min(fl_plot, div_plot);
fill([mach_sorted, fliplr(mach_sorted)], ...
     [ones(1,numel(mach_sorted)), fliplr(fill_y')], ...
     [0.75 1 0.75], 'EdgeColor','none', 'FaceAlpha', 0.35, ...
     'DisplayName', 'Stable region (both mechanisms)');

plot(mach_sorted, fl_plot, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 4, ...
     'DisplayName', sprintf('Flutter margin  (dynamic)%s', ...
                            repmat(' — capped at 50x', fl_capped)));

plot(mach_sorted, div_plot, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 4, ...
     'Color', [0.85 0.15 0.1], ...
     'DisplayName', sprintf('Divergence margin  (static twist)%s', ...
                            repmat(' — capped at 50x', div_capped)));

yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Instability onset  (margin = 1)');

hold off;
grid on; box on;
xlabel({'Mach number', 'M \geq 1.05  (piston theory validity cutoff per Lighthill 1953 / NACA TN 4021)'});
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

%% -----------------------------------------------------------------------
%% 9. Save results
%% -----------------------------------------------------------------------
save(fullfile(figDir, 'flutter.mat'), ...
    'V_fl', 'V_div', 'fl_margin', 'div_margin', ...
    'omega_n', 'f_n', 'Phi_full', 'mesh', ...
    'flightPts', 'fp_crit', 'k_vals', 'flutterResults', ...
    'D_flex', 't', 'rho_lam');

fprintf('\nResults saved to results/flutter.mat\n');
fprintf('Plots saved to results/\n');
fprintf('\n=== Done ===\n');
