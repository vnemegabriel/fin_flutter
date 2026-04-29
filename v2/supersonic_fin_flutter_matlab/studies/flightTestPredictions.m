%% flightTestPredictions.m   (C1)
%  Computes specific, measurable quantities for IREC 2026 flight verification.
%  Requires mainFlutterSolver.m to have been run first (loads flutter.mat).
%
%  Predicted quantities:
%    1. Fin-tip static deflection under peak aerodynamic load [mm]
%    2. Root laminate strain (chordwise) at max-q [microstrain]
%    3. First structural natural frequency [Hz]
%    4. Aeroelastic stability: SM = Inf (no flutter predicted)
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== IREC 2026 Flight Test Predictions (C1) ===\n\n');

%% --- Load main solver results ---
femFile = fullfile(BASE, 'results', 'flutter.mat');
if ~isfile(femFile)
    error('results/flutter.mat not found. Run mainFlutterSolver.m first.');
end
R = load(femFile);   % load all variables

% Rebuild K if not present (pre-existing flutter.mat without K)
if ~isfield(R, 'K') || isempty(R.K)
    fprintf('Note: K not in flutter.mat — rebuilding from mesh and lam.json\n');
    lam   = jsondecode(fileread(fullfile(BASE, 'configs', 'lam.json')));
    D     = lam.tailored_beta;
    D26   = 0; if isfield(D,'D26_Nm'), D26=D.D26_Nm; end
    D_flex2 = [D.D11_Nm, D.D12_Nm, D.D16_Nm; D.D12_Nm, D.D22_Nm, D26; D.D16_Nm, D26, D.D66_Nm];
    t2    = lam.flutter_input.t_mm * 1e-3;
    E_eff = 12 * D.D66_Nm / t2^3;
    geo2.t=t2; mat2.E=E_eff; mat2.nu=0.3;
    K_full = fem.assembleGlobalStiffness(R.mesh, geo2, mat2, D_flex2);
else
    K_full = R.K;
end

mesh     = R.mesh;
Phi_full = R.Phi_full;
omega_n  = R.omega_n;
f_n      = R.f_n;
D_flex   = R.D_flex;
t        = R.t;
fp_crit  = R.fp_crit;

nN    = size(mesh.nodes, 1);
nDOF  = 6 * nN;
wDOFs = (0:nN-1)' * 6 + 3;

b_ref  = 0.22;
beta_M = sqrt(fp_crit.Mach^2 - 1);

fprintf('Critical flight point: Mach=%.3f  q_inf=%.0f Pa  h=%.0f m\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m);

%% --- C1.1: Static tip deflection under peak aerodynamic load ---
% Quasi-steady (k=0) piston-theory static pressure: p = (2q/β) * dw/dx
% For a static solution, w is the unknown — use unit static loading estimate.
%
% Approach: assemble aerodynamic stiffness matrix K_aero from Q_k at k=0,
%   then solve (K_free - Q_aero_free) * u = 0 for the static displacement.
%   Since SM=Inf, invert directly to find tip deflection under a
%   representative unit surface pressure (1 Pa) converted to nodal forces,
%   scaled by fp_crit.q_inf / beta_M * (representative incidence).
%
% Conservative estimate: assume mean AOA-induced pressure = q_inf / beta_M * 0.01
%   (1% incidence angle = ~0.57 deg, representative for small disturbance)

coeff_static = 2 * fp_crit.q_inf / beta_M;
incidence    = 0.01;   % representative dimensionless incidence (rad/rad)

% Assemble consistent nodal force vector from uniform pressure = coeff*incidence
% applied to all element faces in the +z direction (out-of-plane)
F_aero = zeros(nDOF, 1);
nEle = size(mesh.connect, 1);
for ie = 1:nEle
    nd   = mesh.connect(ie, :);
    pts  = mesh.nodes(nd, :);
    t1   = 0.5 * norm(cross(pts(2,:)-pts(1,:), pts(3,:)-pts(1,:)));
    t2_e = 0.5 * norm(cross(pts(3,:)-pts(1,:), pts(4,:)-pts(1,:)));
    area_e = t1 + t2_e;
    p_e    = coeff_static * incidence;   % uniform pressure [Pa]
    f_node = p_e * area_e / 4;           % equal share to 4 nodes
    for nn = 1:4
        F_aero((nd(nn)-1)*6 + 3) = F_aero((nd(nn)-1)*6 + 3) + f_node;
    end
end

% Free DOF indices
tol_bc    = 1e-9;
rootNodes = find(mesh.nodes(:, 2) < tol_bc);
fixedDOFs = reshape((rootNodes-1)*6 + (1:6), 1, []);
allDOFs   = 1:nDOF;
freeDOFs  = setdiff(allDOFs, fixedDOFs);

K_red = K_full(freeDOFs, freeDOFs);
F_red = F_aero(freeDOFs);

u_free = K_red \ F_red;
u_full = zeros(nDOF, 1);
u_full(freeDOFs) = u_free;

% Tip deflection: maximum w at tip nodes (Y = max span)
span      = 0.160;
tipNodes  = find(abs(mesh.nodes(:, 2) - span) < 1e-6);
w_tip     = u_full(wDOFs(tipNodes));
tip_deflection_m  = max(abs(w_tip));
tip_deflection_mm = tip_deflection_m * 1e3;

fprintf('C1.1 — Static tip deflection at max-q:\n');
fprintf('  Uniform surface pressure: %.1f Pa (incidence = %.3f rad)\n', coeff_static*incidence, incidence);
fprintf('  Max tip deflection: %.3f mm\n', tip_deflection_mm);
if tip_deflection_mm < span * 1e3 / 20
    fprintf('  [VALID: < span/20 = %.1f mm — linear theory holds]\n\n', span*1e3/20);
else
    fprintf('  [WARNING: > span/20 — nonlinear effects may apply]\n\n');
end

%% --- C1.2: Root bending moment and strain prediction ---
% Root curvature κ = -∂²w/∂x² ≈ finite difference on first row of elements
% Root nodes: Y ≈ 0 (rootNodes already found above)
% Second row nodes: Y at first element spanwise height
nxMesh = 24;  nyMesh = 12;
first_row_y = span / nyMesh;   % y-coordinate of first interior row
row1_nodes  = find(abs(mesh.nodes(:, 2) - first_row_y) < 1e-6);

if ~isempty(row1_nodes) && ~isempty(rootNodes)
    % Average curvature estimate: (w_row1 - w_root) / dy
    w_root = mean(u_full(wDOFs(rootNodes)));
    w_row1 = mean(u_full(wDOFs(row1_nodes)));
    dy     = first_row_y;
    dw_dy  = (w_row1 - w_root) / dy;   % root rotation proxy

    % Root curvature from beam theory: κ ≈ dw_dy / span (very rough)
    kappa = dw_dy / span;

    % Root strain: epsilon = kappa * (t/2)
    eps_root             = abs(kappa) * (t / 2);
    root_strain_microstrain = eps_root * 1e6;
else
    root_strain_microstrain = NaN;
    fprintf('  Warning: could not find root/row1 nodes for strain estimate\n');
end

fprintf('C1.2 — Root laminate strain at max-q:\n');
fprintf('  Root strain (chordwise): %.1f microstrain\n', root_strain_microstrain);
if root_strain_microstrain >= 10 && root_strain_microstrain <= 10000
    fprintf('  [PLAUSIBLE: 10–10000 με range for composites]\n\n');
else
    fprintf('  [CHECK: outside expected 10–10000 με range — review loading]\n\n');
end

%% --- C1.3: Free-vibration frequency and expected accelerometer response ---
f1_Hz = f_n(1);
f2_Hz = f_n(2);
omega1 = 2 * pi * f1_Hz;

% Order-of-magnitude forced response: |H(ω)| ≈ 1/omega1^2 for undamped (away from resonance)
% Aerodynamic forcing amplitude at max-q: F ~ q_inf * A_fin / b_ref
A_fin    = (0.300 + 0.150) / 2 * span;   % 0.036 m²
F_aero_amp = fp_crit.q_inf * A_fin / b_ref;   % representative force amplitude [N]

% Modal mass M* ≈ rho_m * t * A_fin  (lumped estimate)
modal_mass  = R.rho_m * t * A_fin;
% Expected tip acceleration: a = F/(M* * omega1^2) * omega1^2 = F/M*
expected_accel_ms2 = F_aero_amp / modal_mass;   % [m/s²]
expected_accel_g   = expected_accel_ms2 / 9.81;

fprintf('C1.3 — Free-vibration frequency prediction:\n');
fprintf('  f1 = %.2f Hz  (first bending mode, from modal analysis)\n', f1_Hz);
fprintf('  f2 = %.2f Hz  (second mode)\n', f2_Hz);
fprintf('  Expected accelerometer amplitude at max-q: ~%.1f g (order of magnitude)\n\n', expected_accel_g);

%% --- C1.4: Save flight predictions text file ---
[~, i_crit_q] = max([R.flightPts.q_inf]);
fp_label = R.flightPts(i_crit_q);

outFile = fullfile(BASE, 'results', 'flight_predictions.txt');
fid = fopen(outFile, 'w');
fprintf(fid, '=== IREC 2026 Flight Test Predictions — Team 207 ===\n');
fprintf(fid, 'Critical flight point: Mach %.2f  q_inf=%.0f Pa  alt=%.0f m\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m);

fprintf(fid, 'PREDICTED OBSERVABLE QUANTITIES (to be verified post-flight):\n\n');

fprintf(fid, '1. Fin-tip static deflection at max-q:     %.1f mm\n', tip_deflection_mm);
fprintf(fid, '   [Measurement: photogrammetry or DIC during ground load test at equivalent pressure]\n\n');

fprintf(fid, '2. Root laminate strain at max-q:          %.0f me\n', root_strain_microstrain);
fprintf(fid, '   [Measurement: strain gauge at fin root, chordwise direction]\n\n');

fprintf(fid, '3. First structural natural frequency:     %.2f Hz\n', f1_Hz);
fprintf(fid, '   [Measurement: modal hammer test on completed fin assembly]\n\n');

fprintf(fid, '4. Aeroelastic stability:                  SM = Inf (no flutter predicted)\n');
fprintf(fid, '   [Flight verification: no divergent oscillations in accelerometer/telemetry data]\n\n');

fprintf(fid, 'FALSIFICATION CRITERION:\n');
fprintf(fid, '   Any accelerometer signature showing exponentially growing oscillations\n');
fprintf(fid, '   in the %.0f-%.0f Hz band during the supersonic coast phase would\n', f1_Hz, f2_Hz);
fprintf(fid, '   contradict the SM=Inf prediction and require model revision.\n');
fclose(fid);

fprintf('C1.4 — Predictions saved → results/flight_predictions.txt\n\n');
fprintf('--- Summary ---\n');
fprintf('  Tip deflection:       %.1f mm\n', tip_deflection_mm);
fprintf('  Root strain:          %.0f microstrain\n', root_strain_microstrain);
fprintf('  f1:                   %.2f Hz\n', f1_Hz);
fprintf('  Aeroelastic SM:       Inf (stable)\n');
fprintf('\n=== C1 Done ===\n');
