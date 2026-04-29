%% sensitivityAnalysis.m   (B4)
%  One-at-a-time ±10% material sensitivity + ±2° ply angle scatter.
%  D-matrix recomputed from CLT for each perturbation case.
%  Proves SM=∞ result is robust to ±10% manufacturing variability.
%
%  ar1 layup (from lam.json, 18 physical plies = 4 DB300 biaxial + 5 GA90R per half):
%    Half-stack ordering (outer → mid-plane):
%      Ply 1: DB300 biaxial → sublayers [+45+β, -45+β]
%      Ply 2: GA90R woven   → sublayer  [β]
%      Ply 3: GA90R woven   → sublayer  [β]
%      Ply 4: DB300 biaxial → sublayers [+45+β, -45+β]
%      Ply 5: GA90R woven   → sublayer  [β]
%      Ply 6: DB300 biaxial → sublayers [+45+β, -45+β]
%      Ply 7: GA90R woven   → sublayer  [β]
%      Ply 8: GA90R woven   → sublayer  [β]
%      Ply 9: DB300 biaxial → sublayers [+45+β, -45+β]
%    Bottom half = reverse of top half (symmetric laminate, B=0).
%
%  GA90R woven properties (kc=0.92, Naik & Shembekar 1992):
%    E1_GA = E2_GA = 0.5*(E1+E2)*kc
%    nu12_GA = 2*nu12*E2/(E1+E2)
%    G12_GA  = G12*kc
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== Material Sensitivity Analysis (B4) ===\n\n');

%% --- Nominal material properties (from lam.json) ---
lam   = jsondecode(fileread(fullfile(BASE, 'configs', 'lam.json')));
E1_nom  = lam.material.E1_ud_GPa  * 1e9;
E2_nom  = lam.material.E2_ud_GPa  * 1e9;
G12_nom = lam.material.G12_ud_GPa * 1e9;
nu12_nom = lam.material.nu12_ud;
t       = lam.flutter_input.t_mm * 1e-3;
rho_m   = lam.flutter_input.rho_mat_kgm3;

t_DB = lam.material.t_DB_single_mm * 1e-3;   % single UD sublayer of DB300
t_GA = lam.material.t_GA_mm * 1e-3;          % GA90R woven ply

beta_nom = 20;   % design ply rotation angle [deg]

fprintf('Nominal UD properties:\n');
fprintf('  E1=%.3f GPa  E2=%.3f GPa  G12=%.3f GPa  nu12=%.4f\n', ...
        E1_nom/1e9, E2_nom/1e9, G12_nom/1e9, nu12_nom);
fprintf('  beta=%.0f deg  t_DB=%.4f mm  t_GA=%.4f mm\n\n', ...
        beta_nom, t_DB*1e3, t_GA*1e3);

%% --- B4.1: Perturbation cases (9 total) ---
perturbation_cases = {
    'Nominal (beta=20)',   E1_nom,      E2_nom,      G12_nom,      nu12_nom, beta_nom;
    'E1 x0.90',           E1_nom*0.90, E2_nom,      G12_nom,      nu12_nom, beta_nom;
    'E1 x1.10',           E1_nom*1.10, E2_nom,      G12_nom,      nu12_nom, beta_nom;
    'E2 x0.90',           E1_nom,      E2_nom*0.90, G12_nom,      nu12_nom, beta_nom;
    'E2 x1.10',           E1_nom,      E2_nom*1.10, G12_nom,      nu12_nom, beta_nom;
    'G12 x0.90',          E1_nom,      E2_nom,      G12_nom*0.90, nu12_nom, beta_nom;
    'G12 x1.10',          E1_nom,      E2_nom,      G12_nom*1.10, nu12_nom, beta_nom;
    'beta=18 deg',         E1_nom,      E2_nom,      G12_nom,      nu12_nom, 18;
    'beta=22 deg',         E1_nom,      E2_nom,      G12_nom,      nu12_nom, 22;
};
nCases = size(perturbation_cases, 1);
fprintf('Defined %d perturbation cases.\n\n', nCases);

%% --- Load flight data → critical point ---
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
fprintf('Critical point: Mach=%.3f  q=%.0f Pa\n\n', fp_crit.Mach, fp_crit.q_inf);

%% --- Build mesh + BCs + mass matrix (fixed) ---
cr = 0.300; ct = 0.150; span = 0.160; sweep_deg = 57.4;
nx = 24; ny = 12;
sweep_rad = deg2rad(sweep_deg);

mesh     = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);
nNodes   = size(mesh.nodes, 1);
nDOF_tot = 6 * nNodes;
tol       = 1e-9;
rootNodes = find(mesh.nodes(:, 2) < tol);
fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);

K_dummy          = sparse(nDOF_tot, nDOF_tot);
M_glob           = fem.assembleGlobalMass(mesh, rho_m, t);
[~, M_red, freeDOFs] = fem.applyDirichletBCs(K_dummy, M_glob, fixedDOFs);

nModes = 6;
k_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];
geometry.t = t;
material.nu = 0.3;

%% --- B4.2: Loop over perturbation cases ---
fl_margin_arr    = zeros(nCases, 1);
div_margin_arr   = zeros(nCases, 1);
Q_dyn_mineig_arr = zeros(nCases, 1);
f1_Hz_arr        = zeros(nCases, 1);

fprintf('Running %d perturbation cases...\n', nCases);
fprintf('%-25s  %-8s  %-13s  %-13s  %-15s\n', ...
        'Case', 'f1[Hz]', 'fl_margin', 'div_margin', 'Q_dyn_mineig');
fprintf('%s\n', repmat('-',1,80));

any_finite_fl  = false;
any_finite_div = false;

for ic = 1:nCases
    label  = perturbation_cases{ic, 1};
    E1     = perturbation_cases{ic, 2};
    E2     = perturbation_cases{ic, 3};
    G12    = perturbation_cases{ic, 4};
    nu12   = perturbation_cases{ic, 5};
    beta   = perturbation_cases{ic, 6};

    % Recompute D-matrix from CLT
    D_flex = computeD_ar1(E1, E2, G12, nu12, beta, t_DB, t_GA);

    E_eff      = 12 * D_flex(3,3) / t^3;   % D66-based effective shear modulus
    material.E = E_eff;

    K_glob = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
    [K_red, ~, ~] = fem.applyDirichletBCs(K_glob, M_glob, fixedDOFs);

    [Phi_red, omega_n] = fem.modalAnalysis(K_red, M_red, nModes);
    f_n = omega_n / (2*pi);

    Phi_full = zeros(nDOF_tot, nModes);
    Phi_full(freeDOFs, :) = Phi_red;

    Q_k      = aero.pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                                     fp_crit.a, k_vals, sweep_deg);
    Q_k_norm = Q_k / fp_crit.q_inf;

    [~, ~, res] = stability.solveFlutterPL(omega_n, Q_k_norm, k_vals, fp_crit);

    fl_m   = res.q_flutter_crit / fp_crit.q_inf;
    div_m  = res.q_div_crit    / fp_crit.q_inf;
    mineig = min(res.lam_Qdyn);

    fl_margin_arr(ic)    = fl_m;
    div_margin_arr(ic)   = div_m;
    Q_dyn_mineig_arr(ic) = mineig;
    f1_Hz_arr(ic)        = f_n(1);

    if ~isinf(fl_m),  any_finite_fl  = true; end
    if ~isinf(div_m), any_finite_div = true; end

    if isinf(fl_m),  fl_str  = '    Inf (stable)';  else, fl_str  = sprintf('%13.2f', fl_m);  end
    if isinf(div_m), div_str = '    Inf (stable)';  else, div_str = sprintf('%13.2f', div_m); end

    fprintf('%-25s  %-8.2f  %s  %s  %-15.4e\n', label, f_n(1), fl_str, div_str, mineig);
end
fprintf('%s\n\n', repmat('-',1,80));

if any_finite_fl || any_finite_div
    fprintf('  *** WARNING: one or more perturbed cases yield finite stability margin ***\n');
    fprintf('      Review rows marked with finite margin above.\n\n');
else
    fprintf('  All 9 cases: SM = Inf — design is ROBUST to ±10%% material scatter.\n\n');
end

%% --- B4.3: Plot ---
figDir = fullfile(BASE, 'results');
case_labels = perturbation_cases(:, 1);

figure('Color','w','Visible','off','Position',[50 50 1000 520]);
bar_colors = repmat([0.2 0.5 0.8], nCases, 1);
% Highlight nominal (index 1) in green
bar_colors(1, :) = [0.2 0.7 0.3];
% Highlight any finite margin cases in red
for ic = 1:nCases
    if ~isinf(fl_margin_arr(ic)) || ~isinf(div_margin_arr(ic))
        bar_colors(ic, :) = [0.9 0.2 0.2];
    end
end

b = bar(1:nCases, Q_dyn_mineig_arr, 'FaceColor','flat');
b.CData = bar_colors;
hold on;
yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Zero (neutral)');
hold off;

set(gca, 'XTick', 1:nCases, 'XTickLabel', case_labels, 'XTickLabelRotation', 30);
ylabel('Q_{dyn} minimum eigenvalue  [Pa·(rad/s)^{-2}]');
title({'Material Sensitivity Analysis — Q_{dyn} Minimum Eigenvalue', ...
       'All negative = wash-out stable for all 9 perturbation cases'});

% Annotation
if ~any_finite_fl && ~any_finite_div
    text(0.5, 0.05, 'All SM = \infty  —  Design ROBUST to ±10% scatter', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10, ...
         'Color',[0 0.5 0],'FontWeight','bold');
end
grid on; box on;

saveas(gcf, fullfile(figDir, 'sensitivity.png'));
fprintf('Plot saved → results/sensitivity.png\n');

%% --- B4.4: Text summary ---
Q_min_all = min(Q_dyn_mineig_arr);
Q_max_all = max(Q_dyn_mineig_arr);
all_neg   = all(Q_dyn_mineig_arr < 0);

fid = fopen(fullfile(figDir, 'sensitivity_notes.txt'), 'w');
fprintf(fid, '=== Material Sensitivity Analysis — Robustness Summary ===\n\n');
if any_finite_fl || any_finite_div
    fl_str = 'SOME FINITE VALUES — see table';
else
    fl_str = 'Inf for all 9 cases';
end
fprintf(fid, 'All 9 perturbation cases: SM = %s\n', fl_str);
fprintf(fid, 'Q_dyn_mineig range: [%.4e, %.4e]  (all negative = wash-out robust: %s)\n', ...
        Q_min_all, Q_max_all, mat2str(all_neg));
if ~any_finite_fl && ~any_finite_div
    fprintf(fid, 'Conclusion: Design is ROBUST to ±10%% material scatter.\n');
    fprintf(fid, '  SM = Inf for all 9 perturbation cases.\n');
    fprintf(fid, '  Q_dyn minimum eigenvalue stays negative across all cases,\n');
    fprintf(fid, '  confirming aerodynamic wash-out is not knife-edge.\n');
else
    fprintf(fid, 'Conclusion: Design is SENSITIVE to material scatter.\n');
    fprintf(fid, '  One or more cases yield finite SM — review above table.\n');
end
fclose(fid);
fprintf('Notes saved → results/sensitivity_notes.txt\n');

%% --- Save ---
save(fullfile(figDir, 'sensitivity.mat'), ...
     'perturbation_cases', 'fl_margin_arr', 'div_margin_arr', ...
     'Q_dyn_mineig_arr', 'f1_Hz_arr');
fprintf('Data saved → results/sensitivity.mat\n');
fprintf('\n=== B4 Done ===\n');


%% =========================================================================
%% LOCAL: CLT D-matrix for ar1 laminate
%% =========================================================================
function D_flex = computeD_ar1(E1, E2, G12, nu12, beta_deg, t_DB, t_GA)
% Full CLT bending stiffness matrix for the ar1 layup with global beta rotation.
%
% Layup (from lam.json, 9 physical plies per half-stack):
%   Ply type sequence (outer → mid-plane):
%   DB300, GA90R, GA90R, DB300, GA90R, DB300, GA90R, GA90R, DB300
%   Each DB300 biaxial physical ply → 2 sublayers: [+45+β, -45+β]
%   Each GA90R woven ply → 1 sublayer: [β]  (treated as quasi-isotropic in x-y)
%
% GA90R woven material (kc=0.92, Naik & Shembekar 1992):
%   E1_GA = E2_GA = 0.5*(E1+E2)*0.92
%   nu12_GA = 2*nu12*E2/(E1+E2)
%   G12_GA = G12*0.92

kc     = 0.92;
E1_GA  = 0.5 * (E1 + E2) * kc;
E2_GA  = E1_GA;
nu12_GA = 2.0 * nu12 * E2 / (E1 + E2);
G12_GA = G12 * kc;

% Half-stack sequence: 1=DB300 biaxial, 0=GA90R woven
% Outer → inner (matches lam.json ply order 1..9)
half_type = [1, 0, 0, 1, 0, 1, 0, 0, 1];   % 1=DB300, 0=GA90R

% Build full symmetric sublayer list [angle_deg, t, is_DB]
% Top half (outer → mid) then bottom half (mid → outer = reversed top)
sublayers = build_sublayers(half_type, beta_deg, t_DB, t_GA);
% Mirror: bottom half = reverse of top half
sublayers = [sublayers; flipud(sublayers)];

n_sub   = size(sublayers, 1);
t_total = sum(sublayers(:, 2));

% CLT integration from z = -t/2 to +t/2
z = -t_total / 2;
D = zeros(3, 3);   % only D-matrix needed (symmetric laminate → B=0)

for k = 1:n_sub
    ang_deg = sublayers(k, 1);
    t_k     = sublayers(k, 2);
    is_DB   = sublayers(k, 3);

    if is_DB
        Qbar_k = computeQbar(E1, E2, G12, nu12, ang_deg);
    else
        Qbar_k = computeQbar(E1_GA, E2_GA, G12_GA, nu12_GA, ang_deg);
    end

    z0 = z;
    z1 = z + t_k;
    dz3 = (z1^3 - z0^3) / 3;
    D   = D + Qbar_k * dz3;
    z   = z1;
end

D_flex = D;   % 3×3 bending stiffness [N·m]
end


function subs = build_sublayers(half_type, beta_deg, t_DB, t_GA)
% Convert half-stack ply type sequence to sublayer list [angle, t, is_DB].
% Each DB300 biaxial → 2 sublayers: +45+β and -45+β
% Each GA90R woven   → 1 sublayer:  β
n     = numel(half_type);
subs  = [];
for k = 1:n
    if half_type(k) == 1   % DB300 biaxial
        subs = [subs; +45 + beta_deg, t_DB, 1; ...  %#ok<AGROW>
                      -45 + beta_deg, t_DB, 1];
    else                   % GA90R woven
        subs = [subs; beta_deg, t_GA, 0];            %#ok<AGROW>
    end
end
end


function Qb = computeQbar(E1, E2, G12, nu12, theta_deg)
% Transformed reduced stiffness matrix (3×3 Voigt, in-plane: [11,22,12/66]).
% Qbar maps [eps11, eps22, 2*eps12] → [sig11, sig22, sig12].
nu21 = nu12 * E2 / E1;
denom = 1 - nu12 * nu21;
Q11 = E1  / denom;
Q22 = E2  / denom;
Q12 = nu12 * E2 / denom;
Q66 = G12;

a  = deg2rad(theta_deg);
m  = cos(a);  n = sin(a);
m2 = m^2; n2 = n^2; mn = m*n;

Qb11 = Q11*m2^2 + 2*(Q12 + 2*Q66)*n2*m2 + Q22*n2^2;
Qb22 = Q11*n2^2 + 2*(Q12 + 2*Q66)*n2*m2 + Q22*m2^2;
Qb12 = (Q11 + Q22 - 4*Q66)*n2*m2 + Q12*(m2^2 + n2^2);
Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*n2*m2 + Q66*(n2^2 + m2^2);
Qb16 = (Q11 - Q12 - 2*Q66)*m2*mn - (Q22 - Q12 - 2*Q66)*n2*mn;
Qb26 = (Q11 - Q12 - 2*Q66)*n2*mn - (Q22 - Q12 - 2*Q66)*m2*mn;

% 3×3 format: [11,22,16; 22,22,26; 16,26,66]
Qb = [Qb11, Qb12, Qb16;
      Qb12, Qb22, Qb26;
      Qb16, Qb26, Qb66];
end
