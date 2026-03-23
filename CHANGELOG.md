# CHANGELOG — FalconLAUNCH VI Fin Aeroelastic Scripts

---

## [Unreleased] — 2026-03-21

### `flutterEstimate_v3.py`

#### Fixed — CRITICAL

**F-1 · `ackeret_factor` — wrong AC values replaced by moment-arm formulation**
- **Location:** `flutterEstimate_v3.py:96–128` (function `ackeret_factor`)
- **Before:** `e_sub, e_sup = 0.25, 0.10` — treated AC positions as moment arms; used 0.10c for Ackeret AC (no physical basis).
- **After:** Moment arms computed about the fin elastic axis `x_ea = 0.40c`:
  - `e_sub = x_ea - 0.25` — subsonic AC at 0.25c (thin-airfoil theory)
  - `e_sup = 0.50 - x_ea` — supersonic AC at 0.50c (Ackeret, NACA TM 317)
- **References:** Bisplinghoff, Ashley & Halfman (1955) §5.5; Ackeret, NACA TM 317 (1925)
- **Impact:** Ackeret factor was ~2.55 at design point; corrected to ~1.98 (EA at 0.40c). FSF may change by −22% to −48% depending on EA assumption.

#### Fixed — MAJOR

**F-2 · `compute_flutter` — Ackeret factor now evaluated at self-consistent flutter Mach**
- **Location:** `flutterEstimate_v3.py` function `compute_flutter`, Step 2 block
- **Before:** `f_sup = ackeret_factor(M_rocket)` — evaluated at flight Mach, below the flutter Mach.
- **After:** Fixed-point iteration converges `Mf_iter` such that `Mf = Mf1 * ackeret_factor(Mf)` (≤ 50 iterations, tol = 1e-6). `f_sup` is the factor at the converged flutter Mach.
- **References:** Bisplinghoff, Ashley & Halfman (1955) §5.5

**F-3 · `sweep_factor` — exponent corrected from cos^(−1/2) to 1/cos(Λ)**
- **Location:** `flutterEstimate_v3.py:114–124` (function `sweep_factor`)
- **Before:** `return 1.0 / math.sqrt(math.cos(...))` — factor = cos^(−1/2)(Λ), giving 1.364 at Λ = 57.4°.
- **After:** `return 1.0 / math.cos(...)` — factor = 1/cos(Λ), giving 1.858 at Λ = 57.4°.
- **References:** Bisplinghoff, Ashley & Halfman (1955) §5.5: q_flutter ∝ cos²(Λ) → Vf ∝ cos(Λ).
- **Impact:** Sweep penalty increases from +36% to +86% at design sweep angle.

**F-4 · `flutter_naca4197` — added validation warning for K factor**
- **Location:** `flutterEstimate_v3.py:87–92` (function `flutter_naca4197`)
- **Change:** Added inline warning comment that the linear K formula is unvalidated against NACA TN 4197 Table 2 and ignores taper ratio λ. Validation against Table 2 at (AR, λ) = (0.711, 0.5) is required before relying on FSF results.
- **References:** Martin, NACA TN 4197 (1958) Table 2

#### Fixed — MINOR

**F-5 · `isa` — added troposphere altitude guard**
- **Location:** `flutterEstimate_v3.py:47` (function `isa`)
- **Change:** `assert 0 <= h <= 11000` added as first line. Prevents silent T < 0 errors above 11 km tropopause.

**F-6 · `flutter_naca4197` — documented D16 exclusion**
- **Location:** `flutterEstimate_v3.py` docstring of `flutter_naca4197`
- **Change:** Added note explaining that D16 bend-twist coupling is omitted. Conservative for β = 20° (washout); non-conservative for β > 45° (washin).
- **References:** Weisshaar (1981) J. Aircraft 18(8):669-676

---

### `inplaneG_v5.py`

#### Fixed — MAJOR

**I-1 · `material_props` — crimp knockdown `kc` now applied to woven G12**
- **Location:** `inplaneG_v5.py:103` (function `material_props`)
- **Before:** `'GA': (E1_GA, E1_GA, G12, 0.05, t_GA)` — G12 not degraded by crimp.
- **After:** `'GA': (E1_GA, E1_GA, G12 * kc, nu12_GA, t_GA)` — `kc = 0.92` applied to G12.
- **References:** Naik & Shembekar, J. Compos. Mater. 26(15):2196 (1992)
- **Impact:** Woven G12 reduced by ~8.7%; A66, D66 reduced by ~3–5% for 'ar1' layup.

#### Fixed — MINOR

**I-2 · `material_props` — woven nu12 now derived analytically (replaced hardcoded 0.05)**
- **Location:** `inplaneG_v5.py:102` (function `material_props`)
- **Before:** `0.05` hardcoded with no derivation.
- **After:** `nu12_GA = 2.0 * nu12 * E2 / (E1 + E2)` — CLT cross-ply average for balanced [0/90] stack. Evaluates to ≈ 0.032 at Vf = 0.50.
- **References:** Jones, Mechanics of Composite Materials (1999) §2

**I-3 · `print_layup` — B symmetry check uses dimensionally appropriate reference**
- **Location:** `inplaneG_v5.py:348` (function `print_layup`)
- **Before:** `b_ref = abs(r['A66']) * r['t_total']` — A66·t overestimates reference for thick laminates.
- **After:** `b_ref = abs(r['D66']) / r['t_total']` — D66/t has units [N], same as B matrix, providing a tighter and more physically appropriate normalisation.

---

## Version Notes

These scripts do not carry a standalone semver tag. They are part of the
FalconLAUNCH VI aeroelastic design pipeline at commit state 2026-03-21.
The next tagged release of the full `fin_flutter` package should bump the
minor version to reflect the non-conservative corrections in F-1, F-2, F-3,
and I-1, all of which lower the predicted FSF relative to the previous code.
