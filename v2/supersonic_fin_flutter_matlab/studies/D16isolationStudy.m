%% D16isolationStudy.m   (B2)
%  Isolates the D16 material wash-out contribution by sweeping all beta
%  configurations from lam_sweep.json and comparing against a quasi-isotropic
%  baseline (beta=0 from sweep JSON).
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== D16 Isolation Study (B2) ===\n\n');

%% --- Load lam_sweep.json ---
sweepFile = fullfile(BASE, 'configs', 'lam_sweep.json');
raw     = jsondecode(fileread(sweepFile));
info    = raw.sweep_info;
configs = raw.configs;

t     = info.t_mm * 1e-3;
rho_m = info.rho_mat_kgm3;
nCfg  = numel(configs);

fprintf('Loaded %d beta configs (%.0f to %.0f deg)\n\n', ...
        nCfg, configs(1).beta_deg, configs(end).beta_deg);

%% --- Load flight data → critical flight point ---
opts = detectImportOptions(fullfile(BASE, 'data', 'flight_data.csv'), 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl  = readtable(fullfile(BASE, 'data', 'flight_data.csv'), opts);
h_m  = tbl.altitude_ft * 0.3048;
V    = abs(tbl.Vz_ms);
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
fprintf('Critical point: Mach=%.3f  q=%.0f Pa  h=%.0f m\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m);

%% --- Build mesh + BCs + mass matrix (fixed for all betas) ---
cr = 0.300; ct = 0.150; span = 0.160; sweep_deg = 57.4;
nx = 24; ny = 12;
sweep_rad = deg2rad(sweep_deg);

mesh     = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);
nNodes   = size(mesh.nodes, 1);
nDOF_tot = 6 * nNodes;

tol       = 1e-9;
rootNodes = find(mesh.nodes(:, 2) < tol);
fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);

K_dummy           = sparse(nDOF_tot, nDOF_tot);
M_glob            = fem.assembleGlobalMass(mesh, rho_m, t);
[~, M_red, freeDOFs] = fem.applyDirichletBCs(K_dummy, M_glob, fixedDOFs);

nModes = 6;
k_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];
geometry.t = t;

%% --- B2.1: Sweep all beta configs ---
beta_vec          = zeros(nCfg, 1);
D16_vec           = zeros(nCfg, 1);
fl_margin_vec     = zeros(nCfg, 1);
div_margin_vec    = zeros(nCfg, 1);
Q_dyn_mineig_vec  = zeros(nCfg, 1);
f1_vec            = zeros(nCfg, 1);

fprintf('%-5s  %-9s  %-8s  %-13s  %-13s  %-15s\n', ...
        'beta', 'D16[N.m]', 'f1[Hz]', 'fl_margin', 'div_margin', 'Q_dyn_mineig');
fprintf('%s\n', repmat('-',1,72));

for ic = 1:nCfg
    cfg = configs(ic);
    D26_val = 0;
    if isfield(cfg, 'D26_Nm'), D26_val = cfg.D26_Nm; end
    D_flex = [cfg.D11_Nm, cfg.D12_Nm, cfg.D16_Nm;
              cfg.D12_Nm, cfg.D22_Nm, D26_val;
              cfg.D16_Nm, D26_val,    cfg.D66_Nm];

    E_eff      = 12 * cfg.D66_Nm / t^3;
    material.E  = E_eff;
    material.nu = 0.3;

    K_glob  = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
    [K_red, ~, ~] = fem.applyDirichletBCs(K_glob, M_glob, fixedDOFs);

    [Phi_red, omega_n] = fem.modalAnalysis(K_red, M_red, nModes);
    f_n = omega_n / (2*pi);

    Phi_full = zeros(nDOF_tot, nModes);
    Phi_full(freeDOFs, :) = Phi_red;

    Q_k      = aero.pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                                     fp_crit.a, k_vals, sweep_deg);
    Q_k_norm = Q_k / fp_crit.q_inf;

    [~, ~, res] = stability.solveFlutterPL(omega_n, Q_k_norm, k_vals, fp_crit);

    fl_m     = res.q_flutter_crit / fp_crit.q_inf;
    div_m    = res.q_div_crit    / fp_crit.q_inf;
    mineig   = min(res.lam_Qdyn);

    beta_vec(ic)         = cfg.beta_deg;
    D16_vec(ic)          = cfg.D16_Nm;
    fl_margin_vec(ic)    = fl_m;
    div_margin_vec(ic)   = div_m;
    Q_dyn_mineig_vec(ic) = mineig;
    f1_vec(ic)           = f_n(1);

    if isinf(fl_m),  fl_str  = '    Inf (free)';  else, fl_str  = sprintf('%13.2f', fl_m);  end
    if isinf(div_m), div_str = '    Inf (free)';  else, div_str = sprintf('%13.2f', div_m); end

    fprintf('%-5.0f  %-9.4f  %-8.2f  %s  %s  %-15.4e\n', ...
            cfg.beta_deg, cfg.D16_Nm, f_n(1), fl_str, div_str, mineig);
end
fprintf('%s\n\n', repmat('-',1,72));

%% --- B2.2: Beta=0 case (geometry-only baseline) ---
% Beta=0 is the first config in lam_sweep.json
idx_beta0  = find(beta_vec == 0, 1);
idx_beta20 = find(beta_vec == 20, 1);

if isempty(idx_beta0)
    warning('beta=0 not found in lam_sweep.json; D16 isolation incomplete.');
    Q_dyn_mineig_beta0  = NaN;
    Q_dyn_mineig_beta20 = NaN;
else
    Q_dyn_mineig_beta0  = Q_dyn_mineig_vec(idx_beta0);
    Q_dyn_mineig_beta20 = Q_dyn_mineig_vec(idx_beta20);
end

%% --- B2.3: Tailoring benefit metric ---
fprintf('B2.3 — Tailoring benefit metric:\n');
if ~isnan(Q_dyn_mineig_beta0) && ~isnan(Q_dyn_mineig_beta20)
    tailoring_ratio = Q_dyn_mineig_beta20 / Q_dyn_mineig_beta0;

    fl0  = fl_margin_vec(idx_beta0);
    fl20 = fl_margin_vec(idx_beta20);
    d0   = div_margin_vec(idx_beta0);
    d20  = div_margin_vec(idx_beta20);

    if isinf(fl0),  fl0_str  = 'Inf';  else, fl0_str  = sprintf('%.2f', fl0);  end
    if isinf(fl20), fl20_str = 'Inf';  else, fl20_str = sprintf('%.2f', fl20); end

    fprintf('  beta=0  flutter SM = %s  Q_dyn_mineig = %.4e\n', fl0_str,  Q_dyn_mineig_beta0);
    fprintf('  beta=20 flutter SM = %s  Q_dyn_mineig = %.4e\n', fl20_str, Q_dyn_mineig_beta20);
    fprintf('  Tailoring ratio (beta20/beta0): %.4f\n', tailoring_ratio);

    if tailoring_ratio < 1 && Q_dyn_mineig_beta0 < 0 && Q_dyn_mineig_beta20 < 0
        fprintf('  Interpretation: both geometries are wash-out stable;\n');
        fprintf('    D16 tailoring makes Q_dyn_mineig %.1fx more negative — stronger stabilising force.\n', ...
                abs(tailoring_ratio));
    elseif isinf(fl20) && ~isinf(fl0)
        fprintf('  Interpretation: beta=0 geometry is UNSTABLE (SM=%.2f);\n', fl0);
        fprintf('    D16 tailoring REQUIRED for stability.\n');
    else
        fprintf('  Both configs stable; tailoring ratio = %.4f\n', tailoring_ratio);
    end
else
    tailoring_ratio = NaN;
end

%% --- B2.4: Plot ---
figDir = fullfile(BASE, 'results');
figure('Color','w','Visible','off','Position',[50 50 1000 460]);

fl_plot  = min(fl_margin_vec,  20);
div_plot = min(div_margin_vec, 20);

yyaxis left
plot(beta_vec, fl_plot,  'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b');
hold on;
plot(beta_vec, div_plot, 'b--s','LineWidth', 1.2, 'MarkerFaceColor','b', 'Color',[0.3 0.3 0.9]);
ylabel('Stability margin  (capped at 20)');
set(gca, 'YColor', 'b');

yyaxis right
plot(beta_vec, Q_dyn_mineig_vec, 'r-^', 'LineWidth', 1.5, 'MarkerFaceColor','r');
ylabel('Q_{dyn} minimum eigenvalue  [Pa·(rad/s)^{-2}]');
set(gca, 'YColor', 'r');

if ~isempty(idx_beta20)
    xline(beta_vec(idx_beta20), 'k--', 'LineWidth',1.5, ...
          'Label','\beta=20° (design)', 'LabelVerticalAlignment','top');
end
if ~isempty(idx_beta0)
    xline(beta_vec(idx_beta0), 'g--', 'LineWidth',1.5, ...
          'Label','\beta=0° (geometry-only)', 'LabelVerticalAlignment','bottom');
end

xlabel('\beta [deg]');
title({'D_{16} Isolation: Composite Tailoring vs Geometry-Only Baseline', ...
       sprintf('Critical point: Mach=%.2f  q=%.0f Pa', fp_crit.Mach, fp_crit.q_inf)});
legend({'Flutter margin','Divergence margin','Q_{dyn} mineig'}, ...
       'Location','best','FontSize',8);
grid on; xlim([min(beta_vec)-2, max(beta_vec)+2]);

saveas(gcf, fullfile(figDir, 'D16_isolation.png'));
fprintf('\nPlot saved → results/D16_isolation.png\n');

%% --- Save ---
save(fullfile(figDir, 'D16_isolation.mat'), ...
     'beta_vec', 'D16_vec', 'fl_margin_vec', 'div_margin_vec', ...
     'Q_dyn_mineig_vec', 'f1_vec', 'Q_dyn_mineig_beta0', ...
     'Q_dyn_mineig_beta20', 'tailoring_ratio');
fprintf('Data saved → results/D16_isolation.mat\n');
fprintf('\n=== B2 Done ===\n');
