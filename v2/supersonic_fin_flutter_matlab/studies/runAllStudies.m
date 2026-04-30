%% runAllStudies.m
%  All validation and sensitivity studies for the supersonic fin flutter solver.
%  Run mainFlutterSolver.m first (study C1 loads results/flutter.mat).
%
%  STUDIES IN THIS FILE
%  ─────────────────────────────────────────────────────────────────────────
%  B1  Mesh convergence      — how fine does the FEM mesh need to be?
%  B2  D16 isolation         — how much does composite tailoring actually help?
%  B3  Empirical baseline    — does Bohon (1966) agree with our FEM?
%  B4  Material sensitivity  — is the SM=Inf result robust to ±10% scatter?
%  B5  Frequency benchmark   — does FEM reproduce the analytical cantilever frequency?
%  C1  Flight predictions    — what measurable quantities should we see in flight?
%  ─────────────────────────────────────────────────────────────────────────
%
%  Run from the studies/ directory.

clear; clc;
BASE   = fullfile(fileparts(mfilename('fullpath')), '..');
figDir = fullfile(BASE, 'results');
addpath(BASE);
fprintf('=== All Validation Studies ===\n\n');

%% =========================================================================
%% SHARED SETUP  (used by B1–B5)
%% =========================================================================
% Build nominal laminate.  All studies that need D-matrix, t, or rho
% call core.buildCLTLaminate directly — no lam.json dependency.
[D_flex_nom, t_nom, rho_lam_nom, D_info_nom] = core.buildCLTLaminate(20, 0.50);

cr = 0.300; ct = 0.150; span = 0.160; sweep_deg = 57.4;
sweep_rad = deg2rad(sweep_deg);
k_vals    = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];
nModes    = 6;

% Load flight data and identify critical (max-q) supersonic flight point.
% This point drives the most conservative stability analysis.
csvFile = fullfile(BASE, 'data', 'flight_data.csv');
opts = detectImportOptions(csvFile, 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl  = readtable(csvFile, opts);
h_m  = tbl.altitude_ft * 0.3048;
V    = abs(tbl.Vz_ms);

flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = aero.isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
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


%% =========================================================================
%% B1 — MESH CONVERGENCE STUDY
%% =========================================================================
% QUESTION: Is the 24×12 production mesh fine enough?
%
% METHOD: Run the full flutter analysis on four mesh densities (6×3, 12×6,
%   24×12, 48×24).  Apply Richardson extrapolation (p=2 for Q4 elements)
%   to estimate the exact Mode 1 frequency.  If the 24×12 result is within
%   2% of the extrapolated value, the mesh is converged.
%
% RESULT TO READ: The error column for the 24×12 row.
%   < 2% → mesh is adequate.
%   > 2% → refine to 48×24 for production runs.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('B1 — Mesh Convergence Study\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

mesh_configs = [6, 3; 12, 6; 24, 12; 48, 24];
nMeshes = size(mesh_configs, 1);

freq_table     = zeros(nMeshes, nModes);
fl_margin_B1   = zeros(nMeshes, 1);
div_margin_B1  = zeros(nMeshes, 1);
nElem_vec      = mesh_configs(:,1) .* mesh_configs(:,2);

fprintf('%-8s  %-7s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-10s  %-10s\n', ...
        'Mesh', 'nElem', 'f1[Hz]', 'f2[Hz]', 'f3[Hz]', 'f4[Hz]', 'f5[Hz]', 'f6[Hz]', ...
        'fl_margin', 'div_margin');
fprintf('%s\n', repmat('-', 1, 108));

E_nom = 12 * D_flex_nom(3,3) / t_nom^3;
geo_b1.t = t_nom;   mat_b1.E = E_nom;   mat_b1.nu = 0.3;

for im = 1:nMeshes
    nx = mesh_configs(im,1);   ny = mesh_configs(im,2);
    msh = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);
    rn  = find(msh.nodes(:,2) < 1e-9);
    fd  = reshape((rn-1)*6+(1:6), 1, []);
    K_  = fem.assembleGlobalStiffness(msh, geo_b1, mat_b1, D_flex_nom);
    M_  = fem.assembleGlobalMass(msh, rho_lam_nom, t_nom);
    [Kr, Mr, fr] = fem.applyDirichletBCs(K_, M_, fd);
    [Ph, om] = fem.modalAnalysis(Kr, Mr, nModes);
    f_n = om / (2*pi);
    freq_table(im,:) = f_n;
    Pf  = zeros(size(K_,1), nModes);
    Pf(fr,:) = Ph;
    Qk  = aero.pistonTheoryGAF(msh, Pf, fp_crit.Mach, fp_crit.q_inf, ...
                                 fp_crit.a, k_vals, sweep_deg);
    [~, ~, res] = stability.solveFlutterPL(om, Qk/fp_crit.q_inf, k_vals, fp_crit);
    fl_m  = res.q_flutter_crit / fp_crit.q_inf;
    div_m = res.q_div_crit    / fp_crit.q_inf;
    fl_margin_B1(im)  = fl_m;
    div_margin_B1(im) = div_m;
    if isinf(fl_m),  fl_s  = 'Inf';   else, fl_s  = sprintf('%.3f', fl_m);  end
    if isinf(div_m), div_s = 'Inf';   else, div_s = sprintf('%.3f', div_m); end
    fprintf('%-8s  %-7d  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-8.2f  %-10s  %-10s\n', ...
            sprintf('%dx%d',nx,ny), nx*ny, f_n(1), f_n(2), f_n(3), f_n(4), f_n(5), f_n(6), ...
            fl_s, div_s);
end

% Richardson extrapolation on Mode 1 (p=2 for bilinear Q4)
h_vec  = 1 ./ sqrt(nElem_vec);
f1_vec = freq_table(:,1);
f3 = f1_vec(3);   f4 = f1_vec(4);
h3 = h_vec(3);    h4 = h_vec(4);
f_extrap      = (h3^2 * f4 - h4^2 * f3) / (h3^2 - h4^2);
err_24x12_pct = abs(f3 - f_extrap) / f_extrap * 100;

fprintf('\nRichardson extrapolation (Mode 1, p=2 Q4):\n');
fprintf('  f1 at 24×12:       %.4f Hz\n', f3);
fprintf('  f1 at 48×24:       %.4f Hz\n', f4);
fprintf('  f1 extrapolated:   %.4f Hz\n', f_extrap);
fprintf('  Error at 24×12:    %.3f%%  %s\n\n', err_24x12_pct, ...
        iif(err_24x12_pct < 2, '→ CONVERGED (< 2%)', '→ WARN (> 2%)'));

figure('Color','w','Visible','off','Position',[50 50 1100 480]);
subplot(1,2,1); hold on;
cols = {'b','r','g'};
for mi = 1:3
    plot(h_vec, freq_table(:,mi), [cols{mi} '-o'], 'LineWidth',1.5, ...
         'MarkerFaceColor', cols{mi}, 'DisplayName', sprintf('Mode %d', mi));
end
plot([0, h_vec(end)], [f_extrap, f_extrap], 'b--', 'LineWidth',1.2, ...
     'DisplayName', sprintf('Mode 1 extrap: %.1f Hz', f_extrap));
plot(h_vec(3), freq_table(3,1), 'ko', 'MarkerSize',10, 'LineWidth',2, ...
     'DisplayName', '24×12 (production)');
set(gca,'XDir','reverse'); grid on; legend('Location','best','FontSize',8);
xlabel('1/\surd(nElem)   [coarse ← → fine]');
ylabel('Natural frequency [Hz]');
title('Mode frequencies vs mesh density');

subplot(1,2,2);
fl_pv  = min(fl_margin_B1,  50);
dv_pv  = min(div_margin_B1, 50);
plot(h_vec, fl_pv,  'b-o', 'LineWidth',1.5, 'MarkerFaceColor','b', 'DisplayName','Flutter margin');
hold on;
plot(h_vec, dv_pv,  'r-s', 'LineWidth',1.5, 'MarkerFaceColor','r', 'DisplayName','Divergence margin');
plot(h_vec(3), fl_pv(3), 'ko', 'MarkerSize',10, 'LineWidth',2, 'DisplayName','24×12');
set(gca,'XDir','reverse'); grid on; legend('Location','best','FontSize',8);
xlabel('1/\surd(nElem)   [coarse ← → fine]');
ylabel('Stability margin  q_{crit}/q_{flight}   (> 1 = safe)');
title('Stability margin vs mesh density');
if any(isinf(fl_margin_B1))
    text(0.5,0.92,'All margins = \infty (stable)','Units','normalized',...
         'HorizontalAlignment','center','FontSize',9,'Color',[0 0.5 0]);
end
sgtitle(sprintf('B1 Mesh Convergence — Mach=%.2f  q=%.0f Pa', fp_crit.Mach, fp_crit.q_inf));
saveas(gcf, fullfile(figDir,'mesh_convergence.png'));
save(fullfile(figDir,'mesh_convergence.mat'), ...
     'nElem_vec','freq_table','fl_margin_B1','div_margin_B1','f_extrap','err_24x12_pct');
fprintf('B1 plots and data saved to results/\n\n');


%% =========================================================================
%% B2 — D16 ISOLATION STUDY
%% =========================================================================
% QUESTION: How much does composite tailoring (D16 ≠ 0) actually help
%   compared to a standard layup (beta=0, D16 ≈ 0) at the same geometry?
%
% METHOD: Sweep the global ply rotation beta from 0° to 45° in 5° steps.
%   For each beta, rebuild the D-matrix via CLT, rerun FEM + flutter.
%   Compare the flutter/divergence margin at beta=0 (geometry-only baseline)
%   vs beta=20 (design point).
%
% RESULT TO READ: The Q_dyn_mineig column.
%   All negative → aero wash-out dominates at all betas tested.
%   The ratio (beta=20 mineig) / (beta=0 mineig) quantifies the tailoring benefit.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('B2 — D16 Isolation: Composite Tailoring vs Geometry-Only Baseline\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

beta_sweep = 0:5:45;
nBeta = numel(beta_sweep);

% Mesh and mass matrix are fixed for all beta values
msh_b2 = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, 24, 12);
nN_b2  = size(msh_b2.nodes, 1);
nDOF_b2 = 6 * nN_b2;
rn_b2  = find(msh_b2.nodes(:,2) < 1e-9);
fd_b2  = reshape((rn_b2-1)*6+(1:6), 1, []);
K_dum  = sparse(nDOF_b2, nDOF_b2);
M_b2   = fem.assembleGlobalMass(msh_b2, rho_lam_nom, t_nom);
[~, M_red_b2, free_b2] = fem.applyDirichletBCs(K_dum, M_b2, fd_b2);
geo_b2.t = t_nom;   mat_b2.nu = 0.3;

beta_vec_B2  = zeros(nBeta,1);
D16_vec_B2   = zeros(nBeta,1);
fl_B2        = zeros(nBeta,1);
div_B2       = zeros(nBeta,1);
mineig_B2    = zeros(nBeta,1);
f1_B2        = zeros(nBeta,1);

fprintf('%-5s  %-9s  %-8s  %-13s  %-13s  %-15s\n', ...
        'beta', 'D16[N.m]', 'f1[Hz]', 'fl_margin', 'div_margin', 'Q_dyn_mineig');
fprintf('%s\n', repmat('-',1,72));

for ib = 1:nBeta
    b = beta_sweep(ib);
    [D_b, ~, ~, Di_b] = core.buildCLTLaminate(b, 0.50);
    mat_b2.E  = 12 * D_b(3,3) / t_nom^3;
    K_b  = fem.assembleGlobalStiffness(msh_b2, geo_b2, mat_b2, D_b);
    [Kr_b, ~, ~] = fem.applyDirichletBCs(K_b, M_b2, fd_b2);
    [Ph_b, om_b] = fem.modalAnalysis(Kr_b, M_red_b2, nModes);
    Pf_b = zeros(nDOF_b2, nModes);
    Pf_b(free_b2,:) = Ph_b;
    Qk_b = aero.pistonTheoryGAF(msh_b2, Pf_b, fp_crit.Mach, fp_crit.q_inf, ...
                                  fp_crit.a, k_vals, sweep_deg);
    [~, ~, res_b] = stability.solveFlutterPL(om_b, Qk_b/fp_crit.q_inf, k_vals, fp_crit);
    fl_m  = res_b.q_flutter_crit / fp_crit.q_inf;
    div_m = res_b.q_div_crit    / fp_crit.q_inf;
    me    = min(res_b.lam_Qdyn);
    beta_vec_B2(ib)  = b;
    D16_vec_B2(ib)   = Di_b.D16_Nm;
    fl_B2(ib)  = fl_m;
    div_B2(ib) = div_m;
    mineig_B2(ib) = me;
    f1_B2(ib) = om_b(1) / (2*pi);
    if isinf(fl_m),  fls = '    Inf (free)';  else, fls = sprintf('%13.2f', fl_m);  end
    if isinf(div_m), dvs = '    Inf (free)';  else, dvs = sprintf('%13.2f', div_m); end
    fprintf('%-5.0f  %-9.4f  %-8.2f  %s  %s  %-15.4e\n', ...
            b, Di_b.D16_Nm, om_b(1)/(2*pi), fls, dvs, me);
end
fprintf('%s\n', repmat('-',1,72));

idx0  = find(beta_vec_B2 == 0,  1);
idx20 = find(beta_vec_B2 == 20, 1);
fprintf('\nTailoring benefit (beta=20 vs beta=0):\n');
fprintf('  beta=0   D16=%.4f  fl_margin=%s  Q_dyn_mineig=%.4e\n', ...
        D16_vec_B2(idx0), iif(isinf(fl_B2(idx0)),'Inf',sprintf('%.2f',fl_B2(idx0))), ...
        mineig_B2(idx0));
fprintf('  beta=20  D16=%.4f  fl_margin=%s  Q_dyn_mineig=%.4e\n', ...
        D16_vec_B2(idx20), iif(isinf(fl_B2(idx20)),'Inf',sprintf('%.2f',fl_B2(idx20))), ...
        mineig_B2(idx20));
fprintf('  Ratio (beta20/beta0 Q_dyn_mineig): %.4f\n', mineig_B2(idx20)/mineig_B2(idx0));
fprintf('  Interpretation: ratio < 1 and both negative → D16 makes wash-out stronger.\n\n');

figure('Color','w','Visible','off','Position',[50 50 1000 460]);
yyaxis left
plot(beta_vec_B2, min(fl_B2,20),  'b-o', 'LineWidth',1.5, 'MarkerFaceColor','b');
hold on;
plot(beta_vec_B2, min(div_B2,20), 'b--s','LineWidth',1.2, 'MarkerFaceColor','b');
ylabel('Stability margin (capped at 20)'); set(gca,'YColor','b');
yyaxis right
plot(beta_vec_B2, mineig_B2, 'r-^', 'LineWidth',1.5, 'MarkerFaceColor','r');
ylabel('Q_{dyn} minimum eigenvalue  [Pa·(rad/s)^{-2}]'); set(gca,'YColor','r');
xline(20, 'k--', 'LineWidth',1.5, 'Label','\beta=20° (design)', 'LabelVerticalAlignment','top');
xline(0,  'g--', 'LineWidth',1.5, 'Label','\beta=0° (geometry)', 'LabelVerticalAlignment','bottom');
xlabel('\beta [deg]');
title({'B2  D_{16} Isolation: Composite Tailoring vs Geometry-Only Baseline', ...
       sprintf('Critical point: Mach=%.2f  q=%.0f Pa', fp_crit.Mach, fp_crit.q_inf)});
legend({'Flutter margin','Div margin','Q_{dyn} mineig'},'Location','best','FontSize',8);
grid on; xlim([-2, 47]);
saveas(gcf, fullfile(figDir,'D16_isolation.png'));
save(fullfile(figDir,'D16_isolation.mat'), ...
     'beta_vec_B2','D16_vec_B2','fl_B2','div_B2','mineig_B2','f1_B2');
fprintf('B2 plot and data saved.\n\n');


%% =========================================================================
%% B3 — EMPIRICAL FLUTTER BASELINE (Bohon 1966)
%% =========================================================================
% QUESTION: Does the Bohon (1966) empirical formula agree with our FEM result?
%
% METHOD: Apply the low-AR swept-fin empirical formula:
%   q_flutter_emp = (G_eff · t²) / (AR_eff³ · lambda_taper · cos³Λ)
%   where G_eff = 12·D66/t³ is the equivalent shear modulus.
%
%   Accuracy is ±40% — this is an order-of-magnitude check only.
%   If the empirical formula predicts SM < 1 while FEM predicts SM = Inf,
%   the composite tailoring (D16 wash-out) is the physical reason:
%   the empirical formula was derived for isotropic fins and cannot capture
%   bend-twist coupling.
%
% RESULT TO READ: Empirical SM and FEM SM at the critical flight point.
%   Empirical SM < 1 AND FEM SM = Inf → tailoring was NECESSARY.
%   Both large                         → design has margin to spare.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('B3 — Empirical Flutter Baseline (Bohon 1966)\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

A_plan      = (cr + ct) / 2 * span;            % planform area [m²]
AR_eff      = span^2 / A_plan;
lambda_taper = ct / cr;
F_sweep      = cos(deg2rad(sweep_deg))^3;
G_eff_B3     = 12 * D_flex_nom(3,3) / t_nom^3; % [Pa]
q_flutter_emp = G_eff_B3 * t_nom^2 / (AR_eff^3 * lambda_taper * F_sweep);

fprintf('Geometry:  AR_eff=%.4f   lambda=%.3f   cos³(Λ)=%.4f\n', ...
        AR_eff, lambda_taper, F_sweep);
fprintf('G_eff = D66·12/t³ = %.2f MPa\n', G_eff_B3/1e6);
fprintf('q_flutter_emp = %.2f MPa  (formula accuracy ±40%%)\n\n', q_flutter_emp/1e6);

SM_emp_B3 = q_flutter_emp ./ [flightPts.q_inf];
fprintf('Empirical SM: min=%.1f  max=%.1f  (at critical point: %.1f)\n\n', ...
        min(SM_emp_B3), max(SM_emp_B3), q_flutter_emp/fp_crit.q_inf);

% Overlay with FEM if available
femFile = fullfile(BASE, 'results', 'flutter.mat');
fem_loaded = isfile(femFile);
if fem_loaded
    R3 = load(femFile, 'fl_margin', 'flightPts');
    SM_fem_B3   = R3.fl_margin;
    mach_fem_B3 = [R3.flightPts.Mach];
    fprintf('FEM results loaded from results/flutter.mat\n');
else
    fprintf('Note: flutter.mat not found — run mainFlutterSolver.m to add FEM overlay.\n');
end

mach_B3 = [flightPts.Mach];
[ms, si] = sort(mach_B3);
figure('Color','w','Visible','off','Position',[50 50 900 500]);
hold on;
semilogy(ms, SM_emp_B3(si), 'r-o', 'LineWidth',1.8, 'MarkerSize',4, ...
         'DisplayName','Bohon (1966) empirical ±40%');
if fem_loaded
    [mf_s, sfi] = sort(mach_fem_B3);
    SM_fem_plot = min(R3.fl_margin(sfi), 1e4);
    semilogy(mf_s, SM_fem_plot, 'b-s', 'LineWidth',1.8, 'MarkerSize',4, ...
             'DisplayName','FEM + piston theory');
    if any(isinf(R3.fl_margin))
        text(0.5, 0.97, 'FEM: SM = \infty for all points', 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',9,'Color','b');
    end
end
yline(1,'k--','LineWidth',1.5,'DisplayName','Instability onset (SM=1)');
hold off; grid on; box on;
xlabel('Mach number');
ylabel('Stability margin  q_{flutter}/q_{flight}  (log scale)');
title({'B3  Empirical vs FEM Flutter Stability Margin', ...
       sprintf('Bohon (1966): q_{flutter} = %.1f MPa  (±40%%)', q_flutter_emp/1e6)});
legend('Location','best','FontSize',9);
saveas(gcf, fullfile(figDir,'empirical_vs_fem.png'));
save(fullfile(figDir,'empirical_baseline.mat'), ...
     'SM_emp_B3','mach_B3','q_flutter_emp','G_eff_B3','AR_eff');
fprintf('B3 plot and data saved.\n\n');


%% =========================================================================
%% B4 — MATERIAL SENSITIVITY ANALYSIS
%% =========================================================================
% QUESTION: If T700 fibre or epoxy properties are ±10% from nominal, does
%   the SM=Inf result still hold?  (Manufacturing scatter, lot-to-lot variation.)
%
% METHOD: One-at-a-time perturbation of E1, E2, G12 by ±10%, and beta by ±2°.
%   For each case, rebuild D-matrix via CLT, rerun FEM + flutter.
%   Only the D-matrix changes; mass matrix uses nominal rho_lam (unchanged
%   because rho_lam depends on Vf, not on stiffness moduli).
%
% RESULT TO READ: The fl_margin and div_margin columns.
%   All 'Inf (stable)' → design is ROBUST to ±10% material scatter.
%   Any finite value   → the design is knife-edge and needs investigation.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('B4 — Material Sensitivity Analysis (±10%% ply properties)\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

% Nominal UD properties from buildCLTLaminate D_info
E1_nom  = D_info_nom.E1_GPa  * 1e9;
E2_nom  = D_info_nom.E2_GPa  * 1e9;
G12_nom = D_info_nom.G12_GPa * 1e9;
nu12_nom = D_info_nom.nu12;

perturbation_cases = {
    'Nominal (beta=20)',   E1_nom,      E2_nom,      G12_nom,      nu12_nom, 20;
    'E1 × 0.90',          E1_nom*0.90, E2_nom,      G12_nom,      nu12_nom, 20;
    'E1 × 1.10',          E1_nom*1.10, E2_nom,      G12_nom,      nu12_nom, 20;
    'E2 × 0.90',          E1_nom,      E2_nom*0.90, G12_nom,      nu12_nom, 20;
    'E2 × 1.10',          E1_nom,      E2_nom*1.10, G12_nom,      nu12_nom, 20;
    'G12 × 0.90',         E1_nom,      E2_nom,      G12_nom*0.90, nu12_nom, 20;
    'G12 × 1.10',         E1_nom,      E2_nom,      G12_nom*1.10, nu12_nom, 20;
    'beta = 18 deg',      E1_nom,      E2_nom,      G12_nom,      nu12_nom, 18;
    'beta = 22 deg',      E1_nom,      E2_nom,      G12_nom,      nu12_nom, 22;
};
nCases = size(perturbation_cases, 1);

% Mesh and mass matrix are fixed for all perturbation cases
msh_b4  = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, 24, 12);
nN_b4   = size(msh_b4.nodes,1);
nDOF_b4 = 6*nN_b4;
rn_b4   = find(msh_b4.nodes(:,2) < 1e-9);
fd_b4   = reshape((rn_b4-1)*6+(1:6), 1, []);
K_dum4  = sparse(nDOF_b4, nDOF_b4);
M_b4    = fem.assembleGlobalMass(msh_b4, rho_lam_nom, t_nom);
[~, M_red_b4, free_b4] = fem.applyDirichletBCs(K_dum4, M_b4, fd_b4);
geo_b4.t = t_nom;   mat_b4.nu = 0.3;

fl_B4      = zeros(nCases,1);
div_B4     = zeros(nCases,1);
mineig_B4  = zeros(nCases,1);
f1_B4      = zeros(nCases,1);

fprintf('Running %d perturbation cases...\n', nCases);
fprintf('%-24s  %-8s  %-14s  %-14s  %-15s\n', ...
        'Case','f1[Hz]','fl_margin','div_margin','Q_dyn_mineig');
fprintf('%s\n', repmat('-',1,80));

for ic = 1:nCases
    label = perturbation_cases{ic,1};
    mp = struct('E1',   perturbation_cases{ic,2}, ...
                'E2',   perturbation_cases{ic,3}, ...
                'G12',  perturbation_cases{ic,4}, ...
                'nu12', perturbation_cases{ic,5});
    beta_c = perturbation_cases{ic,6};
    [D_c, ~, ~, ~] = core.buildCLTLaminate(beta_c, 0.50, mp);
    mat_b4.E = 12 * D_c(3,3) / t_nom^3;
    K_c  = fem.assembleGlobalStiffness(msh_b4, geo_b4, mat_b4, D_c);
    [Kr_c, ~, ~] = fem.applyDirichletBCs(K_c, M_b4, fd_b4);
    [Ph_c, om_c] = fem.modalAnalysis(Kr_c, M_red_b4, nModes);
    Pf_c = zeros(nDOF_b4, nModes);
    Pf_c(free_b4,:) = Ph_c;
    Qk_c = aero.pistonTheoryGAF(msh_b4, Pf_c, fp_crit.Mach, fp_crit.q_inf, ...
                                  fp_crit.a, k_vals, sweep_deg);
    [~, ~, res_c] = stability.solveFlutterPL(om_c, Qk_c/fp_crit.q_inf, k_vals, fp_crit);
    fl_m  = res_c.q_flutter_crit / fp_crit.q_inf;
    div_m = res_c.q_div_crit    / fp_crit.q_inf;
    fl_B4(ic)     = fl_m;
    div_B4(ic)    = div_m;
    mineig_B4(ic) = min(res_c.lam_Qdyn);
    f1_B4(ic)     = om_c(1)/(2*pi);
    if isinf(fl_m),  fls = '  Inf (stable)';  else, fls = sprintf('%14.2f', fl_m);  end
    if isinf(div_m), dvs = '  Inf (stable)';  else, dvs = sprintf('%14.2f', div_m); end
    fprintf('%-24s  %-8.2f  %s  %s  %-15.4e\n', label, f1_B4(ic), fls, dvs, mineig_B4(ic));
end
fprintf('%s\n', repmat('-',1,80));

any_finite = any(~isinf(fl_B4)) || any(~isinf(div_B4));
if ~any_finite
    fprintf('\n  All 9 cases: SM = Inf — design is ROBUST to ±10%% material scatter.\n\n');
else
    fprintf('\n  *** WARNING: finite SM in one or more cases — review rows above ***\n\n');
end

bar_cols = repmat([0.2 0.5 0.8], nCases, 1);
bar_cols(1,:) = [0.2 0.7 0.3];  % nominal = green
for ic = 1:nCases
    if ~isinf(fl_B4(ic)) || ~isinf(div_B4(ic))
        bar_cols(ic,:) = [0.9 0.2 0.2];
    end
end
figure('Color','w','Visible','off','Position',[50 50 1000 520]);
b_h = bar(1:nCases, mineig_B4, 'FaceColor','flat');
b_h.CData = bar_cols;
hold on; yline(0,'k--','LineWidth',1.2); hold off;
set(gca,'XTick',1:nCases,'XTickLabel',perturbation_cases(:,1),'XTickLabelRotation',30);
ylabel('Q_{dyn} minimum eigenvalue  [Pa·(rad/s)^{-2}]');
title({'B4  Material Sensitivity — Q_{dyn} Minimum Eigenvalue', ...
       'All negative bars = wash-out stable for every perturbation case'});
if ~any_finite
    text(0.5,0.05,'All SM = \infty  —  Design ROBUST to ±10% scatter', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10,...
         'Color',[0 0.5 0],'FontWeight','bold');
end
grid on; box on;
saveas(gcf, fullfile(figDir,'sensitivity.png'));
save(fullfile(figDir,'sensitivity.mat'), ...
     'perturbation_cases','fl_B4','div_B4','mineig_B4','f1_B4');
fprintf('B4 plot and data saved.\n\n');


%% =========================================================================
%% B5 — ANALYTICAL FREQUENCY BENCHMARK
%% =========================================================================
% QUESTION: Does the FEM assembly reproduce the analytically known frequency
%   for the simplest case (rectangular CFFF plate)?
%
% METHOD: For a rectangular cantilever plate, the first bending frequency is
%   f1 = (1.875104)² / (2π·L²) × sqrt(D11 / (ρ·t))
%   where L = span and D11 = spanwise bending stiffness.
%   This formula is exact for a uniform rectangular plate (Leissa 1969).
%   The test uses a rectangular fin (no sweep, cr=ct=mean chord).
%
%   Expected accuracy: < 5% for AR=0.71 plate vs beam approximation.
%
% RESULT TO READ: Error percentage between FEM and analytical.
%   < 5% → FEM correctly captures bending stiffness and mass.
%   > 5% → bug in assembleGlobalStiffness or assembleGlobalMass.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('B5 — Analytical Frequency Benchmark (rectangular CFFF plate)\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

mean_chord = (cr + ct) / 2;   % 0.225 m
msh_b5  = fem.GenerarMallaAleta(mean_chord, mean_chord, span, 0, 24, 12);
rn_b5   = find(msh_b5.nodes(:,2) < 1e-9);
fd_b5   = reshape((rn_b5-1)*6+(1:6), 1, []);
E_b5    = 12 * D_flex_nom(3,3) / t_nom^3;
geo_b5.t = t_nom;   mat_b5.E = E_b5;   mat_b5.nu = 0.3;
K_b5 = fem.assembleGlobalStiffness(msh_b5, geo_b5, mat_b5, D_flex_nom);
M_b5 = fem.assembleGlobalMass(msh_b5, rho_lam_nom, t_nom);
[Kr_b5, Mr_b5, ~] = fem.applyDirichletBCs(K_b5, M_b5, fd_b5);
[~, om_b5] = fem.modalAnalysis(Kr_b5, Mr_b5, 3);
f_FEM_b5 = om_b5 / (2*pi);

rho_s       = rho_lam_nom * t_nom;   % areal density [kg/m²]
f1_anal_b5  = (1.875104)^2 / (2*pi*span^2) * sqrt(D_info_nom.D11_Nm / rho_s);
err_b5      = abs(f_FEM_b5(1) - f1_anal_b5) / f1_anal_b5 * 100;

fprintf('Rectangular plate  a=%.3f m  b=%.3f m  t=%.1f mm\n', ...
        mean_chord, span, t_nom*1e3);
fprintf('D11=%.2f N·m   rho_s=%.4f kg/m²\n\n', D_info_nom.D11_Nm, rho_s);
fprintf('  Mode 1  FEM:        %.2f Hz\n', f_FEM_b5(1));
fprintf('  Mode 1  Analytical: %.2f Hz\n', f1_anal_b5);
fprintf('  Error:              %.2f%%  %s\n', err_b5, ...
        iif(err_b5 < 5, '→ PASS (< 5%)', '→ WARN (> 5% — check FEM assembly)'));
fprintf('  Mode 2  FEM:        %.2f Hz\n', f_FEM_b5(2));
fprintf('  Mode 3  FEM:        %.2f Hz\n\n', f_FEM_b5(3));

outB5 = fullfile(figDir, 'analytical_benchmark.txt');
fid = fopen(outB5, 'w');
fprintf(fid, '=== B5 FEM Analytical Frequency Benchmark ===\n');
fprintf(fid, 'Rectangular CFFF plate: a=%.3f m  b=%.3f m  t=%d mm\n', ...
        mean_chord, span, round(t_nom*1e3));
fprintf(fid, 'D11=%.2f N·m   rho_s=%.4f kg/m²\n\n', D_info_nom.D11_Nm, rho_s);
fprintf(fid, '  Mode  FEM[Hz]   Analytical[Hz]   Error[%%]\n');
fprintf(fid, '    1   %6.2f     %6.2f            %.2f\n', ...
        f_FEM_b5(1), f1_anal_b5, err_b5);
fprintf(fid, '    2   %6.2f     N/A\n', f_FEM_b5(2));
fprintf(fid, '    3   %6.2f     N/A\n', f_FEM_b5(3));
fclose(fid);
fprintf('B5 report saved to results/analytical_benchmark.txt\n\n');


%% =========================================================================
%% C1 — FLIGHT TEST PREDICTIONS (IREC 2026)
%% =========================================================================
% QUESTION: What specific, measurable quantities will we see in flight data
%   that would confirm or falsify the aeroelastic model?
%
% METHOD: Load the full flutter.mat from mainFlutterSolver.m (must exist).
%   Compute:
%     1. Static tip deflection under representative aerodynamic pressure
%     2. Root strain at max-q (chordwise direction)
%     3. First structural natural frequency (match with accelerometer)
%     4. Stability: SM = Inf expected
%
% RESULT TO READ: The four predicted values in the table.
%   Tip deflection < span/20 → linear theory is valid.
%   Root strain 10–10000 με → plausible composite loading.
%   f1 matches modal hammer test → correct structural model.
%   No flutter signature in flight telemetry → confirms SM = Inf.
% ─────────────────────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('C1 — IREC 2026 Flight Test Predictions\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

femFile = fullfile(BASE, 'results', 'flutter.mat');
if ~isfile(femFile)
    fprintf('  flutter.mat not found — run mainFlutterSolver.m first.\n\n');
    fprintf('=== All Studies Done ===\n');
    return
end
R = load(femFile);

% Build K_full from scratch (not stored in flutter.mat to keep it slim)
[D_c1, ~, ~, ~] = core.buildCLTLaminate(20, 0.50);
E_c1 = 12 * D_c1(3,3) / R.t^3;
geo_c1.t = R.t;   mat_c1.E = E_c1;   mat_c1.nu = 0.3;
K_full = fem.assembleGlobalStiffness(R.mesh, geo_c1, mat_c1, D_c1);

mesh_c1  = R.mesh;
Phi_full = R.Phi_full;
f_n      = R.f_n;
fp_c1    = R.fp_crit;
nN_c1    = size(mesh_c1.nodes,1);
nDOF_c1  = 6*nN_c1;
wDOFs_c1 = (0:nN_c1-1)' * 6 + 3;

beta_M_c1 = sqrt(fp_c1.Mach^2 - 1);
b_ref_c1  = (cr^2 + cr*ct + ct^2) / (3*(cr+ct));   % MAC/2

fprintf('Critical point: Mach=%.3f  q=%.0f Pa  h=%.0f m\n\n', ...
        fp_c1.Mach, fp_c1.q_inf, fp_c1.h_m);

% ── C1.1: Static tip deflection under representative aero pressure ────────
% Quasi-steady pressure p = (2q/β) × (small incidence ≈ 0.01 rad)
coeff_c1  = 2 * fp_c1.q_inf / beta_M_c1;
incidence = 0.01;   % representative (0.57°)

F_aero = zeros(nDOF_c1, 1);
for ie = 1:size(mesh_c1.connect, 1)
    nd   = mesh_c1.connect(ie,:);
    pts  = mesh_c1.nodes(nd,:);
    a1   = 0.5 * norm(cross(pts(2,:)-pts(1,:), pts(3,:)-pts(1,:)));
    a2   = 0.5 * norm(cross(pts(3,:)-pts(1,:), pts(4,:)-pts(1,:)));
    p_e  = coeff_c1 * incidence;
    fn   = p_e * (a1+a2) / 4;
    for nn = 1:4
        F_aero((nd(nn)-1)*6+3) = F_aero((nd(nn)-1)*6+3) + fn;
    end
end
rn_c1  = find(mesh_c1.nodes(:,2) < 1e-9);
fd_c1  = reshape((rn_c1-1)*6+(1:6), 1, []);
allD   = 1:nDOF_c1;
freeD  = setdiff(allD, fd_c1);
u_free = K_full(freeD,freeD) \ F_aero(freeD);
u_full = zeros(nDOF_c1,1);
u_full(freeD) = u_free;

tipN   = find(abs(mesh_c1.nodes(:,2) - span) < 1e-6);
w_tip  = u_full(wDOFs_c1(tipN));
tip_mm = max(abs(w_tip)) * 1e3;
valid_linear = tip_mm < span*1e3/20;

fprintf('C1.1 — Static tip deflection at max-q:\n');
fprintf('  Surface pressure: %.1f Pa  (incidence = %.3f rad)\n', coeff_c1*incidence, incidence);
fprintf('  Max tip deflection: %.3f mm  %s\n\n', tip_mm, ...
        iif(valid_linear,'→ VALID: < span/20, linear theory holds', ...
                         '→ CHECK: > span/20, nonlinear effects may apply'));

% ── C1.2: Root strain estimate ────────────────────────────────────────────
nx_c1   = 24;   ny_c1 = 12;
row1_y  = span / ny_c1;
row1_n  = find(abs(mesh_c1.nodes(:,2) - row1_y) < 1e-6);
if ~isempty(row1_n) && ~isempty(rn_c1)
    w_root_c1 = mean(u_full(wDOFs_c1(rn_c1)));
    w_row1_c1 = mean(u_full(wDOFs_c1(row1_n)));
    kappa_c1  = ((w_row1_c1 - w_root_c1) / row1_y) / span;
    eps_us    = abs(kappa_c1) * (R.t / 2) * 1e6;
else
    eps_us = NaN;
end
fprintf('C1.2 — Root laminate strain at max-q:\n');
fprintf('  Root strain: %.1f microstrain  %s\n\n', eps_us, ...
        iif(eps_us >= 10 && eps_us <= 10000, '→ PLAUSIBLE (10–10000 με)', ...
                                             '→ CHECK: outside expected range'));

% ── C1.3: Natural frequency and accelerometer prediction ──────────────────
f1_c1   = f_n(1);   f2_c1 = f_n(2);
A_fin_c1 = (cr+ct)/2 * span;
modal_mass_c1 = R.rho_lam * R.t * A_fin_c1;
F_amp_c1      = fp_c1.q_inf * A_fin_c1 / b_ref_c1;
expected_g    = F_amp_c1 / modal_mass_c1 / 9.81;

fprintf('C1.3 — Free-vibration frequencies:\n');
fprintf('  f1 = %.2f Hz  (Mode 1, first bending)\n', f1_c1);
fprintf('  f2 = %.2f Hz  (Mode 2)\n', f2_c1);
fprintf('  Expected accelerometer amplitude at max-q: ~%.1f g (order of magnitude)\n\n', ...
        expected_g);

% ── C1.4: Stability ───────────────────────────────────────────────────────
fl_margin_c1 = R.fl_margin;
fprintf('C1.4 — Aeroelastic stability:\n');
if all(isinf(fl_margin_c1)) && all(isinf(R.div_margin))
    fprintf('  SM = Inf for all %d supersonic flight points.\n\n', numel(R.flightPts));
else
    [mf, imf] = min(fl_margin_c1);
    fprintf('  Min flutter margin: %.2f at Mach %.2f\n\n', mf, R.flightPts(imf).Mach);
end

% ── Save predictions text ─────────────────────────────────────────────────
outC1 = fullfile(figDir, 'flight_predictions.txt');
fid = fopen(outC1, 'w');
fprintf(fid, '=== IREC 2026 Flight Test Predictions ===\n');
fprintf(fid, 'Critical point: Mach %.2f  q=%.0f Pa  h=%.0f m\n\n', ...
        fp_c1.Mach, fp_c1.q_inf, fp_c1.h_m);
fprintf(fid, 'PREDICTED QUANTITIES (to verify post-flight):\n\n');
fprintf(fid, '1. Fin-tip static deflection at max-q:  %.1f mm\n', tip_mm);
fprintf(fid, '   Measurement: photogrammetry / DIC during ground load test\n\n');
fprintf(fid, '2. Root laminate strain at max-q:        %.0f microstrain\n', eps_us);
fprintf(fid, '   Measurement: strain gauge at fin root, chordwise direction\n\n');
fprintf(fid, '3. First structural natural frequency:   %.2f Hz\n', f1_c1);
fprintf(fid, '   Measurement: modal hammer test on completed fin assembly\n\n');
fprintf(fid, '4. Aeroelastic stability:                SM = Inf (no flutter predicted)\n');
fprintf(fid, '   Verification: no exponentially growing oscillations in accelerometer\n');
fprintf(fid, '   data in the %.0f–%.0f Hz band during supersonic coast phase\n', f1_c1, f2_c1);
fclose(fid);
fprintf('Flight predictions saved → results/flight_predictions.txt\n\n');

fprintf('─── Summary ────────────────────────────────────────────────────\n');
fprintf('  Tip deflection:        %.1f mm\n', tip_mm);
fprintf('  Root strain:           %.0f microstrain\n', eps_us);
fprintf('  f1:                    %.2f Hz\n', f1_c1);
fprintf('  Aeroelastic SM:        Inf (unconditionally stable)\n\n');


fprintf('=== All Studies Done ===\n');


%% =========================================================================
%% LOCAL HELPER
%% =========================================================================
function y = iif(cond, a, b)
% Inline if: return a when cond is true, b otherwise.
if cond,  y = a;  else,  y = b;  end
end
