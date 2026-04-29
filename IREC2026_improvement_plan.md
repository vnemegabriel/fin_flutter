# IREC 2026 — Solver & Document Improvement Plan
**Team 207 · ITBA · Fin Flutter Aeroelastic Analysis**

> **Purpose**: Task list for Claude Code. Execute blocks in the dependency order shown at the end.
> Each subtask is atomic, self-contained, and includes an acceptance criterion (AC).

---

## Context Block — Read Once, Reuse Always

### Repository layout
```
fin_flutter/v2/supersonic_fin_flutter_matlab/
├── mainFlutterSolver.m          # 8-step orchestrator
├── betaSweepSolver.m            # β-sweep over lam_sweep.json
├── +aero/
│   ├── isaAtmosphere.m
│   └── pistonTheoryGAF.m        # GAF computation (1-pt quadrature)
├── +core/
│   ├── getElementConstitutive.m
│   └── getTransformationMatrix.m
├── +fem/
│   ├── GenerarMallaAleta.m      # mesh builder
│   ├── assembleGlobalStiffness.m
│   ├── assembleGlobalMass.m
│   ├── applyDirichletBCs.m
│   └── modalAnalysis.m
├── +stability/
│   ├── solveFlutterPL.m         # p-L solver (exists but not fully wired)
│   └── loewnerInterpolation.m
├── configs/
│   ├── lam.json                 # flight laminate (β=20°)
│   └── lam_sweep.json           # β-sweep configs
├── data/flight_data.csv
├── results/                     # flutter.mat + PNG outputs
├── studies/                     # CREATE this directory for new scripts
└── tests/
    ├── testFEMAssembly.m
    └── testPLSolver.m
```

### Key physical constants (do NOT re-derive)
```
Fin geometry:   cr=0.300m  ct=0.150m  span=0.160m  Λ=57.4°  t=0.006m
Laminate:       18-ply symmetric  β=20° tailored
D-matrix (N·m): D11=822.01  D22=868.99  D12=264.32  D66=289.76  D16=16.35  D26=11.65
Material:       E1=117.17GPa  E2=7.23GPa  G12=3.47GPa  ν12=0.275  ρ=1580kg/m³
Current mesh:   24×12 Q4 Mindlin-Reissner shell  (nElem=288, nNodes=325)
Natural freqs:  174.62  437.72  986.29  1075.10  1664.02  2051.10  Hz
Modes used:     6
k_vals:         [0, 0.01, 0.05, 0.10, 0.20, 0.50, 1.00]
b_ref:          0.22 m  (reference semi-chord)
Current result: 246 supersonic points  M=1.00–1.85  SM=∞ all
```

### Piston theory validity
```
β = √(M²-1)    → singular at M=1.00 (β=0)
Valid domain:   M ≥ 1.05  (β ≥ 0.312)   per Lighthill (1953), NACA TN 4021
Current bug:    filter uses M_i >= 1.0   → includes the singularity
Fix:            change to M_i >= 1.05
```

### solveFlutterPL.m status
- EXISTS and is mathematically complete (divergence via B̂ eigenvalue, flutter via Q_dyn average)
- NOT called from `mainFlutterSolver.m` — the main file reimplements a subset inline
- Wiring task: replace inline duplication with a call to `solveFlutterPL`

---

## BLOCK A — Critical Solver Fixes

### A1 · Fix M≥1.05 Singularity Filter

**Files**: `mainFlutterSolver.m`, `+aero/pistonTheoryGAF.m`

**Why**: Piston theory has β=0 at M=1.00 → division by zero in `coeff = 2*q_inf/beta`. Currently the flight-point filter passes all M≥1.00, including the singularity. Judging criterion: results near M=1 are scientifically invalid.

#### Subtasks

**A1.1** — `mainFlutterSolver.m`, filter loop (~line 80–100)
- Change: `if M_i >= 1.0` → `if M_i >= 1.05`
- Add comment: `% M >= 1.05 ensures beta = sqrt(M^2-1) >= 0.312; piston theory validity per Lighthill (1953) and NACA TN 4021`
- AC: `sum([flightPts.Mach] >= 1.05)` equals the new supersonic point count; no point has Mach < 1.05.

**A1.2** — `+aero/pistonTheoryGAF.m`, opening guard (~line 20)
- Change existing guard from `if Mach <= 1.0` to:
  ```matlab
  MACH_MIN = 1.05;   % piston theory valid for beta = sqrt(M^2-1) > 0.31
  if Mach < MACH_MIN
      error('pistonTheoryGAF: M=%.4f below validity threshold M=%.2f. beta=%.4f is too small.', ...
            Mach, MACH_MIN, sqrt(max(Mach^2-1, 0)));
  end
  ```
- AC: calling `pistonTheoryGAF` with Mach=1.02 throws an error with the message above.

**A1.3** — `mainFlutterSolver.m`, results reporting block
- Update the `fprintf` summary line to say `M >= 1.05` and report the updated point count.
- Update plot x-axis label: add footnote string `'M \geq 1.05  (piston theory validity cutoff)'`.
- AC: re-run produces updated `flutter_envelope.png` with x-axis starting at ≥1.05.

---

### A2 · Upgrade GAF Quadrature to 2×2 Gauss

**File**: `+aero/pistonTheoryGAF.m`, function `computeModeCoupling`

**Why**: Current implementation uses element centroid (mean of 4 nodal values) for both `w_mean` and `dw/dx`. This is 1-point quadrature — O(h) accuracy. The chordwise gradient `dw/dx` is especially sensitive. 2×2 Gauss is exact for bilinear Q4 shape functions.

**Math reference** (do not re-derive):
```
Bilinear shape functions on [-1,1]²:
  N1=(1-ξ)(1-η)/4   N2=(1+ξ)(1-η)/4   N3=(1+ξ)(1+η)/4   N4=(1-ξ)(1+η)/4

Gauss points (2×2): ξ,η ∈ {-1/√3, +1/√3}, weight = 1 each (total weight = 4)

Jacobian J = [∂x/∂ξ  ∂y/∂ξ; ∂x/∂η  ∂y/∂η],  det(J) gives physical area contribution.

dw/dx in physical space: [dw/dx; dw/dy] = J^{-T} · [∂w/∂ξ; ∂w/∂η]
Only dw/dx (chordwise, x-component) enters the piston-theory pressure.
```

#### Subtasks

**A2.1** — Add a helper sub-function `gaussPoints2x2` at the bottom of `pistonTheoryGAF.m`:
```matlab
function [xi_g, eta_g, w_g] = gaussPoints2x2()
    g = 1/sqrt(3);
    xi_g  = [-g,  g, -g,  g];
    eta_g = [-g, -g,  g,  g];
    w_g   = [ 1,  1,  1,  1];   % weight = 1 per point (sum = 4 for [-1,1]²)
end
```
- AC: `sum(w_g)` == 4; length == 4.

**A2.2** — Rewrite the `computeModeCoupling` inner loop. Replace the current centroid-based approach with:
```
For each element ie:
  pts = mesh.nodes(nd, :)          % [4×3] corner coordinates
  [xi_g, eta_g, w_g] = gaussPoints2x2()
  For each Gauss point g = 1..4:
    N  = [(1-xi_g(g))*(1-eta_g(g)), (1+xi_g(g))*(1-eta_g(g)),
           (1+xi_g(g))*(1+eta_g(g)), (1-xi_g(g))*(1+eta_g(g))] / 4   % [1×4]
    dN_dxi  = [-(1-eta_g(g)),  (1-eta_g(g)),  (1+eta_g(g)), -(1+eta_g(g))] / 4
    dN_deta = [-(1-xi_g(g)), -(1+xi_g(g)),  (1+xi_g(g)),  (1-xi_g(g))] / 4
    J = [dN_dxi; dN_deta] * pts(:,1:2)    % [2×2] Jacobian (x-y plane)
    detJ = det(J)
    dN_dxy = J \ [dN_dxi; dN_deta]        % [2×4] shape fn physical derivatives
    For each mode j:
      w_j_nodes = Phi(wDOFs, j)            % [4×1]
      w_g_val   = N * w_j_nodes            % scalar: w at Gauss point
      dw_dx_g   = dN_dxy(1,:) * w_j_nodes % scalar: dw/dx at Gauss point
      p_j_g = coeff * (1i*k/b_ref * w_g_val + dw_dx_g)
      For each mode i:
        w_i_val = N * Phi(wDOFs, i)
        Qc(i,j) += w_g(g) * w_i_val * p_j_g * detJ
```
- Remove old centroid loop entirely.
- AC: for a uniform planar mode (w=const), `Qc(i,j)` equals `(2*q_inf/beta) * (1i*k/b_ref) * w_const^2 * A_total` to within 0.1%. Run a unit test.

**A2.3** — Add a regression test in `tests/testGAFQuadrature.m`:
```matlab
% Uniform mode (w=1 everywhere): analytical Q_ii = (2q/β)*(ik/b)*A_fin
% A_fin = (cr+ct)/2 * span = (0.3+0.15)/2 * 0.16 = 0.036 m²
% Tolerance: 1%
```
- AC: test passes with < 1% error.

---

### A3 · Wire solveFlutterPL.m into mainFlutterSolver.m

**Files**: `mainFlutterSolver.m` (step 7), `+stability/solveFlutterPL.m`

**Why**: `mainFlutterSolver.m` currently reimplements divergence and flutter logic inline (partially). `solveFlutterPL.m` is the canonical, well-documented implementation. Consolidating eliminates duplication and ensures the p-L method is actually used.

#### Subtasks

**A3.1** — Read `solveFlutterPL.m` signature:
```matlab
[V_fl, V_div, results] = solveFlutterPL(omega_n, Q_k, k_vals, flightConds)
```
`flightConds` is a struct array with fields `.Mach .a .rho .U .q_inf`. Confirm `flightPts` in `mainFlutterSolver.m` already has these fields (it does — check step 2 output).

**A3.2** — In `mainFlutterSolver.m` step 7, replace the inline flutter/divergence block with:
```matlab
%% 7. Flutter + divergence — canonical p-L solver
[V_fl, V_div, flutterResults] = stability.solveFlutterPL( ...
    omega_n, Q_k, k_vals, flightPts);
```
Keep the existing variable names `fl_margin`, `div_margin`, `lam_hat` by reading them from `flutterResults`:
```matlab
fl_margin  = flutterResults.flutter_margin;   % [nF×1]
div_margin = flutterResults.div_margin;       % [nF×1]
lam_hat    = flutterResults.lam_hat;          % [nModes×1] divergence eigenvalues
```
- AC: `mainFlutterSolver.m` runs to completion, produces identical `flutter.mat` fields to the previous version (compare with `isequal` on `fl_margin`, `div_margin` within 1e-10 tolerance).

**A3.3** — Remove duplicate inline flutter/divergence code from `mainFlutterSolver.m` after confirming A3.2 AC passes. Add a comment: `% Flutter/divergence logic lives in +stability/solveFlutterPL.m`.
- AC: `mainFlutterSolver.m` line count decreases; no dead code.

---

## BLOCK B — New Analysis Studies

> All new scripts go in `studies/`. Create the directory if it doesn't exist.
> Each script is standalone: `clear; clc;` at top, loads what it needs, saves results to `results/`.

---

### B1 · Mesh Convergence Study

**New file**: `studies/meshConvergenceStudy.m`

**Why**: Primary validation artefact for the Hoult Award. Proves FEM frequencies and margins converge with mesh refinement; establishes that the 24×12 mesh is well-converged.

**Mesh densities to test**: 6×3, 12×6, 24×12 (current), 48×24 → nElem = 18, 72, 288, 1152

#### Subtasks

**B1.1** — Script structure:
```matlab
clear; clc;
BASE = '..';   % relative to studies/
% Load laminate (lam.json), flight data, build geometry once
% mesh_configs = [6,3; 12,6; 24,12; 48,24]
% For each config:
%   [mesh] = fem.GenerarMallaAleta(nx, nz, ...)
%   [K,M]  = assemble; apply BCs; modal analysis
%   [Q_k]  = aero.pistonTheoryGAF at critical flight point (max q_inf)
%   [~,~,res] = stability.solveFlutterPL(omega_n, Q_k, k_vals, fp_crit)
%   Store: f_n(1:6), fl_margin, div_margin
% Save: results/mesh_convergence.mat
```
- AC: script runs without error; `mesh_convergence.mat` contains `nx_vals`, `nz_vals`, `freq_table` [4×6], `fl_margin_vec` [4×1].

**B1.2** — Richardson extrapolation for Mode 1 frequency:
```
h(i) = 1/sqrt(nElem(i))       % mesh parameter
f_extrap = (h(3)^p * f(4) - h(4)^p * f(3)) / (h(3)^p - h(4)^p)   % p=2 for Q4
error_24x12_pct = abs(f(3) - f_extrap) / f_extrap * 100
```
- AC: Mode 1 error at 24×12 mesh is < 2% vs. extrapolated value.

**B1.3** — Plot: `results/mesh_convergence.png`
- Two subplots: (left) first 3 natural frequencies vs. 1/√nElem with extrapolated dashed line; (right) flutter stability margin vs. 1/√nElem.
- x-axis reversed (coarse left, fine right). Mark 24×12 operating point.
- AC: PNG saved, frequencies show monotone convergence.

---

### B2 · D₁₆ Isolation Parametric Study

**Modify**: `betaSweepSolver.m` (or create `studies/D16isolationStudy.m` if modifying betaSweepSolver risks breaking it)

**Recommended**: create `studies/D16isolationStudy.m` as a clean standalone.

**Why**: Judges noted D₁₆=16.35 N·m is only ~2% of D₁₁. Without isolating the contribution, the "material wash-out" claim is unsubstantiated. This study either proves tailoring adds measurable benefit or honestly reports that geometry dominates.

#### Subtasks

**B2.1** — Script structure:
```matlab
% Load lam_sweep.json (has β = 0..30 in steps)
% For each β config:
%   Build D_flex from JSON; assemble K; modal analysis; compute GAF
%   Run solveFlutterPL at CRITICAL flight point (max q_inf from flight_data.csv)
%   Record: beta_deg, D16_Nm, f1_Hz, fl_margin, div_margin,
%           Q_dyn_mineig (smallest eigenvalue of Q_dyn — the stability indicator)
% Save: results/D16_isolation.mat
```
- AC: `D16_isolation.mat` contains `beta_vec`, `D16_vec`, `fl_margin_vec`, `div_margin_vec`, `Q_dyn_mineig_vec`.

**B2.2** — Explicitly run the β=0° case:
- If `lam_sweep.json` does not include β=0, add a manual D_flex at β=0 using the isotropic approximation: set `D16=0, D26=0`, keep `D11, D22, D12, D66` from the β=0 row or compute them analytically.
- The β=0 case is the "geometry-only" baseline: any finite SM here vs. SM=∞ at β=20° directly quantifies the material wash-out contribution.
- AC: β=0 case executes; `Q_dyn_mineig` at β=0 is stored separately as `Q_dyn_mineig_beta0`.

**B2.3** — Compute the tailoring benefit metric:
```matlab
% If both SM=Inf: compare Q_dyn_mineig values
% More negative mineig at β=20 → stronger wash-out → quantify ratio
tailoring_ratio = Q_dyn_mineig_beta20 / Q_dyn_mineig_beta0;
% If β=0 gives finite SM: SM_beta0 is the direct proof
fprintf('β=0  SM=%.2f  |  β=20°  SM=%.2f  |  Tailoring ratio: %.2f\n', ...)
```
- AC: `tailoring_ratio` printed and saved to `results/D16_isolation.mat`.

**B2.4** — Plot: `results/D16_isolation.png`
- Two y-axes: left = flutter stability margin (cap at 20 for display); right = Q_dyn minimum eigenvalue (always shown).
- Vertical dashed line at β=20° (design point).
- Vertical dashed line at β=0° (geometry-only baseline).
- AC: PNG saved; β=0 and β=20 clearly marked.

---

### B3 · Empirical Flutter Baseline (Bohon/NACA Formula)

**New file**: `studies/empiricalFlutterBaseline.m`

**Why**: Without a baseline, judges cannot assess whether the FEM solver was necessary. If the empirical formula gives a dangerously low margin (finite SM near 1), it proves the high-fidelity approach is needed. If it gives SM=∞ too, it honestly bounds the result.

**Formula** (Bohon 1966, NACA TN 4021, valid for low-AR swept fins):
```
q_flutter_empirical = (G_eff * t²) / (AR_eff³ * λ_taper * F_sweep)

where:
  G_eff    = D66 * 12 / t³       [Pa] — effective shear modulus from CLT
  t        = 0.006 m             — thickness
  AR_eff   = span² / A_planform  — effective aspect ratio
             A_planform = (cr+ct)/2 * span = 0.036 m²
             AR_eff = 0.16² / 0.036 = 0.711
  λ_taper  = ct/cr = 0.5         — taper ratio
  F_sweep  = cos³(Λ) = cos³(57.4°) — sweep factor (approx)

Note: This formula is approximate (±40%). Use it as order-of-magnitude check only.
Document the formula source and limitations explicitly in output.
```

#### Subtasks

**B3.1** — Implement formula, loop over all supersonic flight points:
```matlab
% For each fp in flightPts:
%   q_emp = G_eff * t^2 / (AR_eff^3 * lambda_taper * F_sweep)
%   SM_emp(i) = q_emp / fp.q_inf
% Save: results/empirical_baseline.mat
```
- AC: `SM_emp` vector saved; values are finite (formula always gives a finite q_flutter).

**B3.2** — Overlay plot: `results/empirical_vs_fem.png`
- x-axis: Mach number; y-axis: stability margin
- Line 1: FEM+piston theory SM (from `flutter.mat`)
- Line 2: Bohon empirical SM (from B3.1)
- y-axis log scale if ranges differ greatly
- Legend notes: "FEM (high-fidelity)" and "Bohon (1966) empirical ±40%"
- Horizontal line at SM=1 (instability onset)
- AC: PNG saved; both curves visible; empirical formula gives conservative (lower) SM if the FEM predicts SM=∞.

**B3.3** — Print interpretation text to console and save to `results/empirical_baseline_notes.txt`:
```
Bohon empirical: min SM = X.XX at Mach Y.YY
FEM piston theory: SM = Inf (all points)
Interpretation: [one of the following]
  - "Empirical predicts SM < 1 at some points: FEM analysis was NECESSARY to resolve the design"
  - "Empirical predicts SM = %.1f minimum: FEM confirms and refines with higher confidence"
```
- AC: file saved with correct interpretation string.

---

### B4 · Material Uncertainty / Sensitivity Analysis

**New file**: `studies/sensitivityAnalysis.m`

**Why**: Composite manufacturing has ±5–10% variability in elastic constants. Proves the SM=∞ result is robust to realistic material scatter, not a knife-edge condition.

#### Subtasks

**B4.1** — Define perturbation set (one-at-a-time, ±10%):
```matlab
% Nominal: E1=117.17e9, E2=7.23e9, G12=3.47e9, beta=20 (ply angle tolerance)
% Perturbations (8 cases + 1 nominal = 9 total):
%   E1 ×0.90, E1 ×1.10
%   E2 ×0.90, E2 ×1.10
%   G12×0.90, G12×1.10
%   beta=18°, beta=22°   (±2° ply orientation scatter)
%   Nominal (beta=20°)
```
- AC: `perturbation_cases` cell array with 9 entries defined before any computation.

**B4.2** — For each perturbation case:
- Recompute D-matrix from CLT (re-run the CLT calculation, not just scale D values — since D16 depends nonlinearly on β)
- Assemble K, modal analysis, compute GAF at critical flight point (max q_inf)
- Run `solveFlutterPL`; record `fl_margin`, `div_margin`, `Q_dyn_mineig`, `f1_Hz`
- AC: 9×4 results table printed to console during run.

**B4.3** — Save and plot: `results/sensitivity.mat` + `results/sensitivity.png`
- Bar chart: x-axis = perturbation case label; y-axis = Q_dyn_mineig (most informative when SM=∞ for all)
- If any perturbed case yields finite SM: highlight in red, print warning
- AC: PNG saved; all bars labeled; nominal case highlighted.

**B4.4** — Print robustness conclusion to `results/sensitivity_notes.txt`:
```
All 9 perturbation cases: SM = [Inf/finite values]
Q_dyn_mineig range: [min, max] (all negative = wash-out robust)
Conclusion: Design is [robust / sensitive] to ±10% material scatter.
```
- AC: file saved.

---

### B5 · Analytical Modal Frequency Benchmark

**New file**: `studies/analyticalFrequencyBenchmark.m`

**Why**: Validates FEM assembly correctness without running experiments. Uses the Leissa (1969) analytical solution for a clamped orthotropic rectangular plate.

**Analytical formula** (Leissa 1969, do not re-derive):
```
For clamped rectangular orthotropic plate, a×b (a=chordwise, b=span):
  ω_mn = (β_mn/a²) * sqrt(D11/ρ_s)
where:
  β_mn are tabulated eigenvalue parameters (CLCC boundary: clamped all sides)
  ρ_s = ρ_mat * t   [kg/m²]  — areal density
  For (m=1,n=1): β_11 ≈ 36.0  (for square plate; corrections for aspect ratio)

For simplification: use CFFF (clamped one side, free three sides) plate
  which matches the fin's root-clamped, free-edge BC more accurately.
  CFFF mode (1,1): ω = λ² * sqrt(D11/(ρ_s * a^4))
  λ² ≈ 3.516  (Timoshenko & Woinowsky-Krieger cantilever beam, valid for AR>>1)

Simplified to Euler-Bernoulli cantilever beam (AR = span/mean_chord = 0.16/0.225 = 0.71):
  For AR < 1, plate effects matter. Use 2D Rayleigh-Ritz estimate.
  f1_analytical = (1.875)² / (2π * span²) * sqrt(D11 / (ρ_s))   [Hz]
  where ρ_s = rho_mat * t
```

#### Subtasks

**B5.1** — Simplified rectangular test case geometry:
```matlab
% Use rectangular fin (sweep=0, cr=ct=mean_chord=0.225m, span=0.160m)
% Same t, D_flex, rho_m as flight fin
% Build FEM mesh at 24×12; apply root clamp (same BC)
% Run modal analysis → f_FEM(1:3) [Hz]
```
- AC: `f_FEM` vector of 3 frequencies extracted.

**B5.2** — Analytical cantilever-beam approximation for first mode:
```matlab
rho_s = rho_m * t;           % areal density [kg/m²]
f1_analytical = (1.875104)^2 / (2*pi*span^2) * sqrt(D_flex(1,1) / rho_s);
error_pct = abs(f_FEM(1) - f1_analytical) / f1_analytical * 100;
fprintf('FEM f1 = %.2f Hz | Analytical = %.2f Hz | Error = %.2f%%\n', ...)
```
- AC: error_pct < 5% (beam theory is approximate; 5% is a reasonable tolerance for AR=0.71 plate).

**B5.3** — Save to `results/analytical_benchmark.txt`:
```
=== FEM Analytical Frequency Benchmark ===
Geometry: rectangular CFFF plate  a=0.225m  b=0.160m  t=6mm
Material: D11=822.01 N·m  rho_s=9.48 kg/m²

Mode   FEM [Hz]   Analytical [Hz]   Error [%]
  1    XXX.XX        XXX.XX          X.XX
```
- AC: file saved; all values populated.

---

## BLOCK C — Flight Test Predictions

### C1 · Compute Specific Observable Predictions

**New file**: `studies/flightTestPredictions.m`

**Why**: IREC §2.6.4 explicitly requires "expected results of flight testing." Currently the abstract says only "post-flight structural assessment" — this is not a prediction, it is a description of a process. The solver can compute specific, measurable quantities that could be verified post-flight.

#### Subtasks

**C1.1** — Static tip deflection under peak aerodynamic load:
```matlab
% At critical flight point (max q_inf):
%   Aerodynamic pressure distribution from piston theory (k=0, quasi-steady):
%   p_elem = (2*q_inf/beta) * dw_dx_static   (first-order term only for static)
%   Static solve: K_free * u_static = F_aero
%   where F_aero is assembled from p_elem * area_elem, applied to w-DOFs
%   Tip deflection = max(u_static at tip nodes) in meters
% Report: tip_deflection_mm [mm]
```
- AC: `tip_deflection_mm` computed and printed; must be physically plausible (< span/20 = 8mm for linear theory to hold).

**C1.2** — Root bending moment and strain prediction:
```matlab
% From static solution u_static:
%   Root bending stress: σ = D_flex * κ  where κ = -∂²w/∂x² at root
%   Approximate κ at root: use finite difference of w-DOF values in first row of elements
%   σ_max = D_flex(1,1) * κ * (t/2) / (D_flex(1,1)/E1)  [simplified]
%   Equivalent: ε_root = κ * (t/2)   [microstrain]
% Report: root_strain_microstrain
```
- AC: `root_strain_microstrain` printed; plausible order of magnitude (10–1000 με range for composites).

**C1.3** — Free-vibration frequency prediction (accelerometer-observable):
```matlab
% Already known from modal analysis
% State: "A fin-tip accelerometer triggered during boost should show
%         first bending mode at f1 = 174.62 Hz if attached correctly"
% Compute: expected signal amplitude at max-q using forced response estimate
%   |H(ω)| ≈ 1 / (ω1² - ω²) for undamped single-DOF near resonance
%   At max-q flight, aerodynamic forcing ~ q_inf / b_ref
% Report: f1_Hz, expected_accel_g (order of magnitude)
```
- AC: `f1_Hz` = 174.62 (from modal analysis, just re-state); `expected_accel_g` computed.

**C1.4** — Save all predictions to `results/flight_predictions.txt`:
```
=== IREC 2026 Flight Test Predictions — Team 207 ===
Critical flight point: Mach X.XX  q_inf=XXXX Pa  alt=XXXXX m

PREDICTED OBSERVABLE QUANTITIES (to be verified post-flight):
1. Fin-tip static deflection at max-q:     XX.X mm
   [Measurement: photogrammetry or DIC during ground load test at equivalent pressure]

2. Root laminate strain at max-q:          XXX με
   [Measurement: strain gauge at fin root, chordwise direction]

3. First structural natural frequency:     174.62 Hz
   [Measurement: modal hammer test on completed fin assembly]

4. Aeroelastic stability:                  SM = Inf (no flutter predicted)
   [Flight verification: no divergent oscillations in accelerometer/telemetry data]

FALSIFICATION CRITERION:
   Any accelerometer signature showing exponentially growing oscillations
   in the 174–440 Hz band during the supersonic coast phase would
   contradict the SM=Inf prediction and require model revision.
```
- AC: file saved with all four predictions populated with numerical values.

---

## BLOCK D — Document Revisions

> Execute ONLY after Block A (fixes) and at least B1, B2, B5 (key validation results) are complete.
> Use actual numerical values from the study outputs — do not invent numbers.

---

### D1 · Rewrite Extended Abstract

**File**: `Team207_ITBA_Extended_Abstract_2026IREC.docx`

**Constraint**: ≥500 words, ≤2 pages, AIAA style, NO equations/figures/tables.

**Note on template**: The IREC rules require the official ESRA AIAA Word template from esrarocket.org. If the template can be downloaded, use it. If not, maintain the current formatting which follows AIAA conventions.

#### Subtasks

**D1.1** — Add mission context paragraph (insert at end of §I Introduction):
```
"The vehicle operates with a [motor class] motor. The simulated trajectory reaches
a peak dynamic pressure of [q_max] Pa at Mach [M_at_qmax] and [altitude] m altitude.
The supersonic coast phase lasts approximately [duration] seconds, during which the
fins must remain aeroelastically stable. Piston theory is applied for Mach numbers
at or above 1.05, corresponding to the validity domain β = √(M²−1) ≥ 0.31 per
Lighthill (1953) and NACA TN 4021, yielding [updated_count] analyzed flight points."
```
→ Replace bracketed values with actual numbers from the solver output.

**D1.2** — Revise §III Theoretical Framework, piston theory paragraph:
- Add: "...validity requiring Mach ≥ 1.05. The analysis excludes the transonic region near Mach 1.0 where the Prandtl-Glauert factor β approaches zero and the linearized theory breaks down..."
- This directly addresses the singularity critique.

**D1.3** — Add validation subsection to §IV Results (rename from just "Results" to "Results and Validation"):
- Add sentences (not a section header — no subheadings in AIAA extended abstract):
  "A mesh convergence study across four refinement levels from 18 to 1152 elements confirms that the 24×12 mesh used for all analyses achieves a Mode 1 frequency error below [X]% relative to Richardson-extrapolated values. The cantilever-beam analytical solution for Mode 1 yields [f_analytical] Hz, against the FEM result of [f_FEM] Hz — a [Y]% difference consistent with the plate-versus-beam geometry correction."

**D1.4** — Add D₁₆ isolation result:
- In §V Physical Interpretation, after describing the dual mechanism, add:
  "To isolate the material wash-out contribution, the analysis was repeated with a quasi-isotropic laminate (β=0°, D₁₆≈0) retaining the identical fin geometry. [Result: either 'The stability margin decreased from infinity to X, directly quantifying the tailoring benefit' OR 'Both configurations retain infinite stability margin; however, the minimum eigenvalue of the dynamic aerodynamic stiffness matrix is [ratio]× more negative for the tailored design, indicating a proportionally larger stabilizing force.'"]"

**D1.5** — Rewrite §VI Conclusions, add specific flight test predictions:
- Replace the current vague "post-flight structural assessment" with:
  "Flight testing at IREC 2026 will be used to verify three specific predictions derived from the computational model: (1) first bending natural frequency of [f1] Hz, verifiable by pre-flight modal hammer test; (2) fin-tip deflection of [X] mm under a representative ground aerodynamic load test; and (3) absence of any exponentially growing oscillations in the [f1]–[f2] Hz band in accelerometer telemetry during the supersonic coast phase. Any sustained, growing oscillation in this band would constitute evidence of flutter onset and require model revision."

**D1.6** — Review and tighten all undefined terms:
- First use of "aeroelastic tailoring": add "(the engineering of laminate fiber orientation to generate desired coupling between structural deformation modes)"
- First use of "wash-out": add "(a twist response that reduces aerodynamic angle of attack under bending load)"
- Remove or define "positive-definite" if retained in the text.
- AC: reread the full document. Every technical term is defined on first use.

---

### D2 · Update Presentation Slides

**File**: `Team207_ITBA_Presentation_2026IREC.pptx` (regenerate via `make_slides.js`)

#### Subtasks

**D2.1** — Update slide 9 (Results) with corrected point count (after M≥1.05 filter):
- Change "246" KPI card to the new count.
- Update Mach range if lower bound shifts from 1.00 to 1.05.

**D2.2** — Add slide 9b (new slide after Results): "Model Validation"
- Left panel: mesh convergence — table of 4 mesh sizes with Mode 1 frequency and % change
- Right panel: analytical benchmark — FEM vs. cantilever beam, error percentage
- Slide title: "Validation: Mesh Convergence & Analytical Benchmark"
- Position: between current slides 9 and 10; renumber subsequent slides.

**D2.3** — Add slide 9c (new slide): "D₁₆ Isolation: Quantifying the Material Wash-Out"
- Left panel: β sweep plot showing SM (or Q_dyn_mineig) vs. β with β=0 and β=20° marked
- Right panel: text explanation with the tailoring ratio or SM comparison
- Slide title: "Composite Tailoring: Isolating the D₁₆ Contribution"

**D2.4** — Update slide 3 (Problem Statement) — add empirical baseline context:
- Add one bullet or text sentence: "The Bohon (1966) empirical formula predicts a minimum SM of [X] — the FEM analysis [confirms / refines this / contradicts this]."

**D2.5** — Add "Flight Test Predictions" content to slide 12 (Conclusions) or as a dedicated slide 12b:
- Three-item list with specific values from C1.4
- Frame as: "Falsifiable predictions — not assumptions"

**D2.6** — Fix slide 5 (Theoretical Framework) piston theory pillar text:
- Add: "Valid for M ≥ 1.05  (β = √(M²−1) ≥ 0.31)"
- Acknowledges the domain limitation proactively.

**D2.7** — Update total slide count in footer (`n / 13` → `n / 15` after adding 2 new slides).
- AC: all slides regenerated, QA pass (render to images via LibreOffice + pdftoppm, inspect all 15 slides for overflow/overlap).

---

## Execution Order & Dependencies

```
A1 ──┐
A2 ──┼──► Run mainFlutterSolver.m (clean baseline) ──► B5 ──► B1 ──► B2 ──► B3 ──► B4 ──► C1
A3 ──┘                                                                                    │
                                                                                          ▼
                                                                                    D1 ──► D2
```

### Priority tiers

| Tier | Blocks | Rationale |
|---|---|---|
| 🔴 Critical (do first) | A1, A2, A3 | Fix correctness bugs; all subsequent results depend on a correct solver |
| 🟡 High (before documents) | B5, B1, B2 | Primary validation; D₁₆ isolation; needed for D1 text |
| 🟢 Medium | B3, B4, C1 | Strengthen argument; flight predictions |
| 🔵 Last | D1, D2 | Documents use actual numbers from studies |

---

## Acceptance Criteria Summary

| Task | Pass condition |
|---|---|
| A1 | No flight point with Mach < 1.05 in `flightPts`; `pistonTheoryGAF(1.02,...)` throws error |
| A2 | Uniform-mode unit test passes with < 1% error; `testGAFQuadrature.m` green |
| A3 | `mainFlutterSolver.m` output matches previous `flutter.mat` within 1e-10; line count reduced |
| B1 | `mesh_convergence.mat` saved; Mode 1 convergence error < 2% at 24×12 mesh |
| B2 | β=0 case runs; `tailoring_ratio` or SM comparison printed; `D16_isolation.png` saved |
| B3 | `empirical_vs_fem.png` saved; interpretation string correct |
| B4 | All 9 perturbation cases run; `sensitivity.png` saved; robustness note saved |
| B5 | Analytical benchmark error < 5%; `analytical_benchmark.txt` saved |
| C1 | `flight_predictions.txt` saved with 4 populated predictions and falsification criterion |
| D1 | All 6 subtasks complete; undefined terms eliminated; word count 500–1100 words; ≤ 2 pages |
| D2 | 15 slides rendered; all overflow/overlap checks pass; new slides visible and readable |
