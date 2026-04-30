%% runAllTests.m
%  Full validation suite for the supersonic fin flutter solver.
%  Runs 9 tests covering atmosphere, FEM, aerodynamics, and stability.
%
%  HOW TO READ THIS FILE
%  Each test block states:
%    WHAT   — the quantity being tested
%    WHY    — why a failure here is a physical problem, not just a number
%    PASS criterion — the specific threshold
%
%  Run from the supersonic_fin_flutter_matlab/ directory.

clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
fprintf('=== Flutter Solver Validation Suite ===\n\n');
PASS = true;

%% ── Shared setup ────────────────────────────────────────────────────────────
% Build laminate via embedded CLT (no lam.json dependency).
[D_flex, t, rho_lam, D_info] = core.buildCLTLaminate(20, 0.50);
E_eff  = 12 * D_flex(3,3) / t^3;
geometry.t  = t;
material.E  = E_eff;
material.nu = 0.3;

% Coarse 6×3 mesh — fast for tests; real solver uses 24×12
mesh      = fem.GenerarMallaAleta(0.300, 0.150, 0.160, deg2rad(57.4), 6, 3);
rootNodes = find(mesh.nodes(:,2) < 1e-9);
fixedDOFs = reshape((rootNodes-1)*6 + (1:6), 1, []);

K = fem.assembleGlobalStiffness(mesh, geometry, material, D_flex);
M = fem.assembleGlobalMass(mesh, rho_lam, t);
[K_red, M_red, freeDOFs] = fem.applyDirichletBCs(K, M, fixedDOFs);

nModes = 4;
[Phi_red, omega_test] = fem.modalAnalysis(K_red, M_red, nModes);
Phi_full = zeros(size(K,1), nModes);
Phi_full(freeDOFs,:) = Phi_red;


%% =========================================================================
%% T1 & T2 — ISA atmosphere
%% =========================================================================
% WHAT  : The International Standard Atmosphere model (ISA 1975).
%         This function is called at every flight point to convert altitude
%         to density, speed of sound, temperature and pressure.
%
% WHY   : If air density or speed of sound are wrong, every Mach number,
%         dynamic pressure, and flutter margin is wrong.  The checks use
%         values tabulated in ISO 2533:1975 to within rounding.
%
% PASS  : Within 0.1% of ISO 2533 tabulated values.

fprintf('─── Atmosphere (T1–T2) ─────────────────────────────────────────────\n');
fprintf('Tests the ISA model used to convert altitude → density, Mach, q_inf.\n\n');

[rho0, a0, T0, P0] = aero.isaAtmosphere(0);
check(abs(rho0 - 1.225)  < 0.001,  'T1a rho  at h=0 m',   rho0,  1.225,  PASS); PASS = ans;
check(abs(a0   - 340.29) < 0.50,   'T1b a    at h=0 m',   a0,    340.29, PASS); PASS = ans;
check(abs(T0   - 288.15) < 0.01,   'T1c T    at h=0 m',   T0,    288.15, PASS); PASS = ans;
check(abs(P0   - 101325) < 1,      'T1d P    at h=0 m',   P0,    101325, PASS); PASS = ans;

[~, a11, T11, ~] = aero.isaAtmosphere(11000);
check(abs(T11 - 216.65) < 0.01,    'T2a T    at h=11km',  T11,   216.65, PASS); PASS = ans;
check(abs(a11 - 295.07) < 0.50,    'T2b a    at h=11km',  a11,   295.07, PASS); PASS = ans;

h_vec = [0; 5000; 11000; 15000];
[rho_v, ~, ~, ~] = aero.isaAtmosphere(h_vec);
check(length(rho_v) == 4,           'T2c vectorised length',length(rho_v), 4, PASS); PASS = ans;

fprintf('\n');


%% =========================================================================
%% T3 — Stiffness matrix positive-definite after BCs
%% =========================================================================
% WHAT  : The reduced stiffness matrix K_red (free DOFs only, after
%         removing the clamped root) must have all positive eigenvalues.
%
% WHY   : A non-positive-definite K means the structure has a zero-energy
%         (mechanism) mode — some rigid-body motion or hourglass mode has
%         not been eliminated.  Flutter eigenvalues computed from such a
%         system are meaningless.  The tolerance allows for numerical noise
%         at the ~1e-8 level relative to the largest eigenvalue.
%
% PASS  : min eigenvalue > −1e-8 × max eigenvalue.

fprintf('─── FEM Stiffness (T3) ─────────────────────────────────────────────\n');
fprintf('Checks K_red has no zero-energy modes after root clamping.\n\n');

ev_K   = real(eig(full(K_red)));
tol_pd = 1e-8 * max(abs(ev_K));
check(min(ev_K) > -tol_pd, 'T3 K_red positive-definite', min(ev_K), -tol_pd, PASS);
PASS = ans;
fprintf('\n');


%% =========================================================================
%% T4 — Mass matrix positive-definite
%% =========================================================================
% WHAT  : All eigenvalues of M_red must be strictly positive.
%
% WHY   : A zero mass eigenvalue means a DOF has no inertia — the modal
%         analysis will produce a spurious infinite-frequency mode or divide
%         by zero inside the eigensolver.  This catches bugs in the lumped
%         mass assembly (e.g. missing rotational inertia terms).
%
% PASS  : All eigenvalues > 0 (exact — lumped mass has no numerical noise).

fprintf('─── FEM Mass (T4) ──────────────────────────────────────────────────\n');
fprintf('Checks M_red has no zero-mass DOFs (would produce spurious modes).\n\n');

eig_M = real(eig(full(M_red)));
check(all(eig_M > 0), 'T4 M_red positive-definite', min(eig_M), 0, PASS);
PASS = ans;
fprintf('\n');


%% =========================================================================
%% T5 — First natural frequency within 30% of beam-strip estimate
%% =========================================================================
% WHAT  : The first structural natural frequency from FEM vs the
%         Euler-Bernoulli cantilever-strip formula.
%
%         Beam-strip formula:  f1 = (1.875)² / (2π·L²) × sqrt(D22 / (ρ·t))
%
% WHY   : This formula over-predicts the FEM result by 15–25% because it
%         ignores (a) free leading/trailing-edge boundaries, (b) taper, and
%         (c) sweep coupling.  Nevertheless, a >30% discrepancy flags gross
%         errors in D22, ρ, t, or span — the four most flutter-critical
%         parameters.  The 30% tolerance is deliberately wide to remain
%         valid across mesh densities.
%
% PASS  : |f1_FEM − f1_beam| / f1_beam < 30%.

fprintf('─── Natural Frequency (T5) ─────────────────────────────────────────\n');
fprintf('Compares FEM Mode 1 to Euler-Bernoulli cantilever-strip formula.\n\n');

f1_fem  = omega_test(1) / (2*pi);
span    = 0.160;
f1_anal = (1.875104)^2 / (2*pi*span^2) * sqrt(D_info.D22_Nm / (rho_lam * t));
err_pct = abs(f1_fem - f1_anal) / f1_anal * 100;
check(err_pct < 30, 'T5 f1 vs beam-strip  (< 30%)', err_pct, 30, PASS);
PASS = ans;
fprintf('     FEM = %.2f Hz   beam-strip ≈ %.2f Hz\n\n', f1_fem, f1_anal);


%% =========================================================================
%% T6 — Mode shapes mass-orthonormal
%% =========================================================================
% WHAT  : The modal matrix Phi_red satisfies  Phi' M Phi = I  (identity)
%         to numerical precision.
%
% WHY   : The p-L flutter solver works in modal coordinates; it assumes the
%         modal mass matrix is the identity.  If orthonormality fails,
%         the generalised stiffness and aerodynamic force matrices are
%         incorrectly scaled, and flutter speeds are wrong.
%
% PASS  : Frobenius norm of (Phi'·M·Phi − I) < 1e-6.

fprintf('─── Mode Shape Orthonormality (T6) ─────────────────────────────────\n');
fprintf('Verifies Phi'' * M * Phi = I (required by modal-coordinate flutter solver).\n\n');

ortho_err = norm(Phi_red' * M_red * Phi_red - eye(nModes), 'fro');
check(ortho_err < 1e-6, 'T6 ||Phi''MPhi − I||_F', ortho_err, 1e-6, PASS);
PASS = ans;
fprintf('\n');


%% =========================================================================
%% T7 — GAF matrix Hermitian
%% =========================================================================
% WHAT  : Q_k(:,:,ki) must satisfy  Q = Q'  (Hermitian) for every
%         reduced frequency ki.
%
% WHY   : The flutter eigenvalue problem K_ae = K_modal − q·Q_dyn must have
%         real eigenvalues (so that stability margins are real numbers).
%         This requires Q_dyn to be Hermitian.  pistonTheoryGAF enforces
%         Q = (Q + Q')/2 explicitly — this test verifies that symmetrisation
%         worked correctly.
%
% PASS  : max|Q(i,j) − conj(Q(j,i))| < 1e-10.

fprintf('─── GAF Hermitian (T7) ─────────────────────────────────────────────\n');
fprintf('Checks the aerodynamic force matrix Q_k is Hermitian at all frequencies.\n');
fprintf('(Hermitian Q_k ensures the flutter eigenvalue problem has real solutions.)\n\n');

[rho_a, a_a, ~, ~] = aero.isaAtmosphere(4000);
Mach_test = 2.0;
q_test    = 0.5 * rho_a * (Mach_test * a_a)^2;
k_vals    = [0, 0.01, 0.1, 0.5, 1.0];

Q_k = aero.pistonTheoryGAF(mesh, Phi_full, Mach_test, q_test, a_a, k_vals, 57.4);

herm_errs = zeros(1, length(k_vals));
for ki = 1:length(k_vals)
    herm_errs(ki) = max(max(abs(Q_k(:,:,ki) - Q_k(:,:,ki)')));
end
check(max(herm_errs) < 1e-10, 'T7 Q_k Hermitian (all k)', max(herm_errs), 1e-10, PASS);
PASS = ans;
fprintf('\n');


%% =========================================================================
%% T8 — GAF imaginary/real ratio at k→0
%% =========================================================================
% WHAT  : At k=0 (purely quasi-steady), the imaginary part of Q_k should
%         be exactly zero.  At k=0.01 it should be much smaller than the
%         real part.
%
% WHY   : The imaginary part of Q_k is the aerodynamic damping: it is
%         proportional to k (reduced frequency).  At k=0 there is no
%         oscillation, so there can be no damping.  A large Im/Re ratio at
%         low k would mean spurious aerodynamic damping is contaminating the
%         quasi-steady divergence analysis.
%
% PASS  : ||Im(Q_k(:,:,2))||_F / ||Re(Q_k(:,:,2))||_F < 0.5  at k=0.01.

fprintf('─── GAF Im/Re ratio (T8) ────────────────────────────────────────────\n');
fprintf('Checks that imaginary (damping) part of Q_k is small at low frequency.\n\n');

Qr    = real(Q_k(:,:,2));   % k=0.01
Qi    = imag(Q_k(:,:,2));
ratio = norm(Qi,'fro') / (norm(Qr,'fro') + eps);
check(ratio < 0.5, 'T8 Im/Re at k=0.01  (< 0.5)', ratio, 0.50, PASS);
PASS = ans;
fprintf('\n');


%% =========================================================================
%% T9 — Flutter/divergence solver end-to-end
%% =========================================================================
% WHAT  : The p-L solver must run without error and return finite or Inf
%         values for V_fl and V_div that are physically plausible.
%
% WHY   : This is an integration test.  It catches interface bugs between
%         GAF, modal analysis, and the eigenvalue tracker — e.g. wrong
%         normalisation of Q_k, wrong k=0 divergence branch, or complex
%         arithmetic errors that break the eigenvalue sweep.
%
%         For the AR1 laminate at Mach 2, we expect either:
%           (a) V_fl = Inf / V_div = Inf   — wash-out design, both stable
%           (b) V_fl > 0, V_div > 0       — finite speeds, must be positive
%         In either case, V_fl must not be NaN and must be > 0 if finite.
%
% PASS  : ~isnan(V_fl)  AND  (isinf(V_fl) OR V_fl > 0)
%         ~isnan(V_div) AND  (isinf(V_div) OR V_div > 0)

fprintf('─── Flutter solver end-to-end (T9) ─────────────────────────────────\n');
fprintf('Runs the full p-L eigenvalue solver and checks outputs are valid.\n\n');

fp_t.Mach  = Mach_test;   fp_t.a   = a_a;      fp_t.rho   = rho_a;
fp_t.U     = Mach_test * a_a;
fp_t.q_inf = q_test;

Q_k_norm = Q_k / q_test;
[V_fl_t, V_div_t, ~] = stability.solveFlutterPL(omega_test, Q_k_norm, k_vals, fp_t);

fl_ok  = ~isnan(V_fl_t)  && (isinf(V_fl_t)  || V_fl_t  > 0);
div_ok = ~isnan(V_div_t) && (isinf(V_div_t) || V_div_t > 0);

check(fl_ok,  'T9a V_flutter valid  (>0 or Inf, not NaN)', V_fl_t,  0, PASS);
PASS = ans;
check(div_ok, 'T9b V_diverge valid  (>0 or Inf, not NaN)', V_div_t, 0, PASS);
PASS = ans;

if isinf(V_fl_t)
    fprintf('     Flutter speed:   STABLE (Inf) — expected for wash-out AR1 laminate\n');
else
    fprintf('     Flutter speed:   %.1f m/s  (Mach %.2f)\n', V_fl_t, V_fl_t/a_a);
end
if isinf(V_div_t)
    fprintf('     Divergence speed: STABLE (Inf)\n\n');
else
    fprintf('     Divergence speed: %.1f m/s\n\n', V_div_t);
end


%% =========================================================================
%% Summary
%% =========================================================================
fprintf('════════════════════════════════════════════════════════════════════\n');
if PASS
    fprintf('  All tests PASSED.\n');
else
    fprintf('  One or more tests FAILED — review output above.\n');
end
fprintf('════════════════════════════════════════════════════════════════════\n');


%% =========================================================================
%% Helper
%% =========================================================================
function result = check(cond, name, val, expected, current_pass)
if isnumeric(val) && isscalar(val)
    if isinf(val),  valStr = 'Inf';
    elseif isnan(val), valStr = 'NaN';
    else,           valStr = sprintf('%.4g', val);
    end
else
    valStr = num2str(val);
end
if ischar(expected) || isstring(expected)
    expStr = expected;
elseif isnumeric(expected) && isscalar(expected)
    expStr = sprintf('%.4g', expected);
else
    expStr = num2str(expected);
end
if cond
    fprintf('  %-44s  PASS  (got %s)\n', name, valStr);
else
    fprintf('  %-44s  FAIL  (got %s, expected %s)\n', name, valStr, expStr);
end
result = current_pass && cond;
end
