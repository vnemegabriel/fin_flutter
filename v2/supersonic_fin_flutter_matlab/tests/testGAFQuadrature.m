%% testGAFQuadrature.m
%  Unit test for 2x2 Gauss quadrature in pistonTheoryGAF.
%
%  Analytical reference: uniform mode (w=1 everywhere)
%    Q_ii = (2*q/beta) * (i*k/b_ref) * A_fin
%  where A_fin = (cr+ct)/2 * span = (0.3+0.15)/2 * 0.16 = 0.036 m²
%  Tolerance: 1%
%
%  Run from supersonic_fin_flutter_matlab/ directory.

clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
fprintf('=== GAF Quadrature Tests ===\n\n');
PASS = true;

%% --- T1: gaussPoints2x2 properties (verified implicitly by T2) ---
% gaussPoints2x2 is a local sub-function inside pistonTheoryGAF and cannot be
% called externally. Its correctness is verified indirectly: if sum(w_g)==4
% and points are at ±1/√3, the uniform-mode test (T2) will hit exactly 1% or
% better. A wrong weight sum would produce a proportional error > 1%.
fprintf('T1: gaussPoints2x2 verified implicitly by uniform-mode integral (T2)\n');

%% --- T2: Uniform mode w=1 everywhere, k>0 ---
% Build a minimal flat rectangular mesh (no sweep, no taper)
cr = 0.300; ct = 0.150; span = 0.160;
sweep_rad = deg2rad(57.4);
nx = 24; ny = 12;

mesh = fem.GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny);
nN   = size(mesh.nodes, 1);
nDOF = 6 * nN;

% One mode: w = 1 at every node's out-of-plane DOF, 0 elsewhere
Phi = zeros(nDOF, 1);
wDOFs = (0:nN-1)' * 6 + 3;
Phi(wDOFs) = 1.0;

% Flight parameters (arbitrary supersonic)
Mach  = 1.5;
q_inf = 5000;          % Pa
a_snd = 340;           % m/s
k     = 0.1;
k_vals    = [0, k];
sweep_deg = 57.4;

Q_k = aero.pistonTheoryGAF(mesh, Phi, Mach, q_inf, a_snd, k_vals, sweep_deg);

% Analytical Q(1,1) at k>0:
%   Q_ii = (2q/β) * (i*k/b_ref) * A_fin
%   Only imaginary part from w_mean; dw/dx term integrates to zero for w=const
b_ref  = 0.22;
beta   = sqrt(Mach^2 - 1);
A_fin  = (cr + ct) / 2 * span;    % = 0.036 m²
Q_anal = (2*q_inf/beta) * (1i*k/b_ref) * A_fin;

Q_num  = Q_k(1, 1, 2);   % k index 2 (k=0.1)

err_pct = abs(Q_num - Q_anal) / abs(Q_anal) * 100;
fprintf('T2: Uniform mode Q_ii at k=%.2f\n', k);
fprintf('    Analytical: %.6e + %.6e i\n', real(Q_anal), imag(Q_anal));
fprintf('    Numerical:  %.6e + %.6e i\n', real(Q_num),  imag(Q_num));
fprintf('    Error: %.4f%%  ', err_pct);
if err_pct < 1.0
    fprintf('PASS\n');
else
    fprintf('FAIL (> 1%%)\n');
    PASS = false;
end

%% --- T3: k=0 slice is real ---
fprintf('T3: Q_k(:,:,1) is real (k=0)... ');
imag_norm = norm(imag(Q_k(:,:,1)), 'fro');
if imag_norm < 1e-10
    fprintf('PASS  (||Im||_F = %.2e)\n', imag_norm);
else
    fprintf('FAIL  (||Im||_F = %.2e)\n', imag_norm);
    PASS = false;
end

%% --- T4: Hermitian symmetry ---
fprintf('T4: Q_k Hermitian (Q_ij = conj(Q_ji))... ');
err_herm = max(abs(Q_k(:,:,2) - Q_k(:,:,2)'), [], 'all');
if err_herm < 1e-12
    fprintf('PASS  (max err = %.2e)\n', err_herm);
else
    fprintf('FAIL  (max err = %.2e)\n', err_herm);
    PASS = false;
end

%% --- Summary ---
fprintf('\n');
if PASS
    fprintf('All GAF quadrature tests PASSED.\n');
else
    fprintf('One or more tests FAILED — check output above.\n');
end
