%% analyticalFrequencyBenchmark.m   (B5)
%  Validates FEM assembly by comparing Mode 1 frequency to the Euler-Bernoulli
%  cantilever-beam analytical solution for a rectangular CFFF plate.
%
%  Geometry: rectangular fin  a=0.225 m (mean chord)  b=0.160 m (span)
%  Same thickness, D-matrix, and density as the flight fin.
%  BC: root clamp (Y=0), all other edges free.
%
%  Run from the studies/ directory.

clear; clc;
BASE = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(BASE);
fprintf('=== Analytical Frequency Benchmark (B5) ===\n\n');

%% --- Load laminate properties ---
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

fprintf('D11=%.2f  D22=%.2f  D66=%.2f  N·m\n', D.D11_Nm, D.D22_Nm, D.D66_Nm);
fprintf('t=%.1f mm   rho=%.0f kg/m³\n\n', t*1e3, rho_m);

%% --- B5.1: Rectangular fin mesh (sweep=0, cr=ct=mean_chord=0.225 m) ---
mean_chord = 0.225;   % (cr+ct)/2 = (0.300+0.150)/2
span       = 0.160;
nx         = 24;
ny         = 12;

mesh = fem.GenerarMallaAleta(mean_chord, mean_chord, span, 0, nx, ny);

tol       = 1e-9;
rootNodes = find(mesh.nodes(:, 2) < tol);
fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);

K = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
M = fem.assembleGlobalMass(mesh, rho_m, t);
[K_red, M_red, ~] = fem.applyDirichletBCs(K, M, fixedDOFs);

nModes = 3;
[~, omega_n] = fem.modalAnalysis(K_red, M_red, nModes);
f_FEM = omega_n / (2*pi);

fprintf('FEM natural frequencies (rectangular CFFF plate):\n');
fprintf('  Mode 1: %.2f Hz\n', f_FEM(1));
fprintf('  Mode 2: %.2f Hz\n', f_FEM(2));
fprintf('  Mode 3: %.2f Hz\n\n', f_FEM(3));

%% --- B5.2: Analytical cantilever-beam approximation (Mode 1) ---
% Euler-Bernoulli cantilever:  f1 = (1.875104)^2 / (2*pi*L^2) * sqrt(D11/rho_s)
% Using D11 (chordwise bending stiffness) is appropriate when the plate is
% bent in the span direction by the first cantilever mode.
% For AR=0.71 the plate correction vs beam is < 10%; 5% tolerance used.
rho_s       = rho_m * t;                    % areal density [kg/m²]
f1_analytical = (1.875104)^2 / (2*pi*span^2) * sqrt(D.D11_Nm / rho_s);

error_pct = abs(f_FEM(1) - f1_analytical) / f1_analytical * 100;
fprintf('B5.2 — Analytical benchmark:\n');
fprintf('  rho_s   = %.4f kg/m²\n', rho_s);
fprintf('  f1_anal = %.2f Hz   (Euler-Bernoulli cantilever, L=%.3f m)\n', f1_analytical, span);
fprintf('  f1_FEM  = %.2f Hz\n', f_FEM(1));
fprintf('  Error   = %.2f%%\n\n', error_pct);

if error_pct < 5
    fprintf('  PASS (< 5%% tolerance for AR=%.2f plate vs beam)\n\n', ...
            span / mean_chord);
else
    fprintf('  WARN (> 5%% — check mesh or D-matrix)\n\n');
end

%% --- B5.3: Save report ---
outFile = fullfile(BASE, 'results', 'analytical_benchmark.txt');
fid = fopen(outFile, 'w');
fprintf(fid, '=== FEM Analytical Frequency Benchmark ===\n');
fprintf(fid, 'Geometry: rectangular CFFF plate  a=%.3fm  b=%.3fm  t=%dmm\n', ...
        mean_chord, span, round(t*1e3));
fprintf(fid, 'Material: D11=%.2f N·m  rho_s=%.2f kg/m²\n', D.D11_Nm, rho_s);
fprintf(fid, '\n');
fprintf(fid, 'Mode   FEM [Hz]   Analytical [Hz]   Error [%%]\n');
fprintf(fid, '  1    %6.2f        %6.2f          %5.2f\n', ...
        f_FEM(1), f1_analytical, error_pct);
fprintf(fid, '  2    %6.2f        N/A               N/A\n', f_FEM(2));
fprintf(fid, '  3    %6.2f        N/A               N/A\n', f_FEM(3));
fprintf(fid, '\n');
fprintf(fid, 'Note: Analytical formula valid for Mode 1 only (Euler-Bernoulli cantilever).\n');
fprintf(fid, 'AR = %.2f — plate effects cause FEM to deviate from beam theory by up to 10%%.\n', ...
        span / mean_chord);
fclose(fid);
fprintf('Report saved to results/analytical_benchmark.txt\n');
fprintf('\n=== B5 Done ===\n');
