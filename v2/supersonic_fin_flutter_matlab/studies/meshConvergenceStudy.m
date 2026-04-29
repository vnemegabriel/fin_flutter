%% meshConvergenceStudy.m   (B1)
%  Mesh convergence study: natural frequencies and flutter margin vs mesh density.
%  Mesh densities: 6×3, 12×6, 24×12 (current), 48×24  → nElem = 18, 72, 288, 1152
%  Richardson extrapolation on Mode 1 frequency (p=2 for Q4 elements).
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== Mesh Convergence Study (B1) ===\n\n');

%% --- Load laminate and geometry ---
lam   = jsondecode(fileread(fullfile(BASE, 'configs', 'lam.json')));
D     = lam.tailored_beta;
D26   = 0;
if isfield(D, 'D26_Nm'), D26 = D.D26_Nm; end
D_flex = [D.D11_Nm, D.D12_Nm, D.D16_Nm;
          D.D12_Nm, D.D22_Nm, D26;
          D.D16_Nm, D26,      D.D66_Nm];

t     = lam.flutter_input.t_mm * 1e-3;
rho_m = lam.flutter_input.rho_mat_kgm3;
E_eff = 12 * D.D66_Nm / t^3;

geometry.t  = t;
material.E  = E_eff;
material.nu = 0.3;

cr        = 0.300; ct  = 0.150; span = 0.160; sweep_deg = 57.4;
sweep_rad = deg2rad(sweep_deg);

%% --- Load flight data → critical flight point ---
opts = detectImportOptions(fullfile(BASE, 'data', 'flight_data.csv'), 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl  = readtable(fullfile(BASE, 'data', 'flight_data.csv'), opts);

h_m = tbl.altitude_ft * 0.3048;
V   = abs(tbl.Vz_ms);
flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = aero.isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
    % M >= 1.05 ensures beta = sqrt(M^2-1) >= 0.312; piston theory validity per Lighthill (1953) and NACA TN 4021
    if M_i >= 1.05
        fp.Mach=M_i; fp.a=a_i; fp.rho=rho_i;
        fp.U=V(i); fp.q_inf=0.5*rho_i*V(i)^2;
        fp.h_m=h_m(i); fp.time=tbl.time_s(i);
        if isempty(flightPts), flightPts=fp;
        else,                  flightPts(end+1)=fp; %#ok<AGROW>
        end
    end
end
[~, i_crit] = max([flightPts.q_inf]);
fp_crit = flightPts(i_crit);
fprintf('Critical flight point: Mach=%.3f  q=%.0f Pa  h=%.0f m\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m);

%% --- B1.1: Mesh sweep ---
mesh_configs = [6, 3; 12, 6; 24, 12; 48, 24];
nMeshes = size(mesh_configs, 1);
nModes  = 6;
k_vals  = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];

freq_table    = zeros(nMeshes, nModes);
fl_margin_vec = zeros(nMeshes, 1);
div_margin_vec = zeros(nMeshes, 1);
nx_vals       = mesh_configs(:, 1);
nz_vals       = mesh_configs(:, 2);
nElem_vec     = nx_vals .* nz_vals;

fprintf('%-8s  %-7s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-12s  %-12s\n', ...
        'Mesh', 'nElem', 'f1[Hz]', 'f2[Hz]', 'f3[Hz]', 'f4[Hz]', 'f5[Hz]', 'f6[Hz]', ...
        'fl_margin', 'div_margin');
fprintf('%s\n', repmat('-', 1, 108));

for im = 1:nMeshes
    nx = mesh_configs(im, 1);
    ny = mesh_configs(im, 2);

    mesh = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);

    tol       = 1e-9;
    rootNodes = find(mesh.nodes(:, 2) < tol);
    fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);

    K = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
    M_mat = fem.assembleGlobalMass(mesh, rho_m, t);
    [K_red, M_red, freeDOFs] = fem.applyDirichletBCs(K, M_mat, fixedDOFs);

    [Phi_red, omega_n] = fem.modalAnalysis(K_red, M_red, nModes);
    f_n = omega_n / (2*pi);
    freq_table(im, :) = f_n;

    Phi_full = zeros(size(K, 1), nModes);
    Phi_full(freeDOFs, :) = Phi_red;

    Q_k      = aero.pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                                     fp_crit.a, k_vals, sweep_deg);
    Q_k_norm = Q_k / fp_crit.q_inf;

    [~, ~, res] = stability.solveFlutterPL(omega_n, Q_k_norm, k_vals, fp_crit);
    fl_m  = res.q_flutter_crit / fp_crit.q_inf;
    div_m = res.q_div_crit / fp_crit.q_inf;

    fl_margin_vec(im)  = fl_m;
    div_margin_vec(im) = div_m;

    if isinf(fl_m),  fl_str  = 'Inf';  else, fl_str  = sprintf('%.3f', fl_m);  end
    if isinf(div_m), div_str = 'Inf';  else, div_str = sprintf('%.3f', div_m); end

    fprintf('%-8s  %-7d  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-12s  %-12s\n', ...
            sprintf('%dx%d', nx, ny), nx*ny, f_n(1), f_n(2), f_n(3), f_n(4), f_n(5), f_n(6), ...
            fl_str, div_str);
end

%% --- B1.2: Richardson extrapolation on Mode 1 ---
% h(i) = 1/sqrt(nElem(i)), p=2 for Q4 elements
h_vec    = 1 ./ sqrt(nElem_vec);
f1_vec   = freq_table(:, 1);
p_order  = 2;

% Use the two finest meshes (index 3 and 4: 24×12 and 48×24)
h3 = h_vec(3); h4 = h_vec(4);
f3 = f1_vec(3); f4 = f1_vec(4);
f_extrap = (h3^p_order * f4 - h4^p_order * f3) / (h3^p_order - h4^p_order);
error_24x12_pct = abs(f3 - f_extrap) / f_extrap * 100;

fprintf('\nB1.2 — Richardson extrapolation (Mode 1, p=2 Q4):\n');
fprintf('  f1 at 24x12:   %.4f Hz\n', f3);
fprintf('  f1 at 48x24:   %.4f Hz\n', f4);
fprintf('  f_extrapolated: %.4f Hz\n', f_extrap);
fprintf('  Error at 24x12: %.3f%%\n\n', error_24x12_pct);

if error_24x12_pct < 2
    fprintf('  PASS (< 2%% convergence criterion)\n\n');
else
    fprintf('  WARN (> 2%% — mesh may not be converged)\n\n');
end

%% --- B1.3: Plot ---
figDir = fullfile(BASE, 'results');

h_plot   = 1 ./ sqrt(nElem_vec);
h_extrap = 0;   % extrapolated (h=0, infinite mesh)

figure('Color','w','Visible','off','Position',[50 50 1100 480]);

% Left: first 3 natural frequencies vs 1/sqrt(nElem)
subplot(1,2,1);
colors = {'b','r','g'};
hold on;
for mi = 1:3
    plot(h_plot, freq_table(:, mi), [colors{mi} '-o'], 'LineWidth', 1.5, ...
         'MarkerFaceColor', colors{mi}, 'DisplayName', sprintf('Mode %d', mi));
end
% Richardson extrapolated dashed line for Mode 1
plot([0, h_plot(end)], [f_extrap, f_extrap], 'b--', 'LineWidth', 1.2, ...
     'DisplayName', sprintf('Mode 1 extrapolated: %.1f Hz', f_extrap));
% Mark 24×12 operating point
plot(h_plot(3), freq_table(3,1), 'ko', 'MarkerSize', 10, 'LineWidth', 2, ...
     'DisplayName', '24×12 (operating)');
set(gca, 'XDir', 'reverse');
xlabel('1/\surd(nElem)   [coarse ← → fine]');
ylabel('Natural frequency [Hz]');
title('Mode frequencies vs mesh density');
legend('Location','best','FontSize',8); grid on;

% Right: flutter stability margin vs 1/sqrt(nElem)
subplot(1,2,2);
fl_plot_vec  = min(fl_margin_vec,  50);
div_plot_vec = min(div_margin_vec, 50);
plot(h_plot, fl_plot_vec,  'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b', ...
     'DisplayName','Flutter margin');
hold on;
plot(h_plot, div_plot_vec, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor','r', ...
     'DisplayName','Divergence margin');
plot(h_plot(3), fl_plot_vec(3), 'ko', 'MarkerSize', 10, 'LineWidth', 2, ...
     'DisplayName','24×12 (operating)');
set(gca, 'XDir', 'reverse');
xlabel('1/\surd(nElem)   [coarse ← → fine]');
ylabel('Stability margin  q_{crit}/q_{flight}   (> 1 = safe)');
title('Aeroelastic stability margin vs mesh density');
legend('Location','best','FontSize',8); grid on;
if any(isinf(fl_margin_vec))
    text(0.5, 0.92, 'All margins = \infty (stable)', 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',9,'Color',[0 0.5 0]);
end

sgtitle(sprintf('Mesh Convergence Study  —  Critical point: Mach=%.2f  q=%.0f Pa', ...
        fp_crit.Mach, fp_crit.q_inf), 'FontSize', 12);

saveas(gcf, fullfile(figDir, 'mesh_convergence.png'));
fprintf('Plot saved → results/mesh_convergence.png\n');

%% --- Save ---
save(fullfile(figDir, 'mesh_convergence.mat'), ...
     'nx_vals', 'nz_vals', 'nElem_vec', 'freq_table', ...
     'fl_margin_vec', 'div_margin_vec', 'f_extrap', 'error_24x12_pct');
fprintf('Data saved → results/mesh_convergence.mat\n');
fprintf('\n=== B1 Done ===\n');
