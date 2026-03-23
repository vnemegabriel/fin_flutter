# ERRORS.md — Code Review Report: `flutterEstimate_v3.py` & `inplaneG_v5.py`
## FalconLAUNCH VI Supersonic Fin Aeroelastic Design

**Date:** 2026-03-18
**Reviewer:** Claude Code (automated static analysis + theory check)
**Files reviewed:** `last_files/flutterEstimate_v3.py`, `last_files/inplaneG_v5.py`, `last_files/laminate.json`
**Theory references:** NACA TN 4197 (Martin 1958), Ackeret (NACA TM 317), Bisplinghoff/Ashley/Halfman (1955), Jones (1999), Halpin-Tsai (1969), Weisshaar (1981)

---

## Executive Summary

Both scripts are structurally sound and mathematically coherent for subsonic-regime structural analysis (CLPT, ABD matrix, ISA atmosphere). However, **three errors in `flutterEstimate_v3.py` directly affect the supersonic flutter speed prediction** and together produce an **over-optimistic (non-conservative) flutter boundary**, which is a safety concern for a supersonic rocket fin. One inconsistency exists in `inplaneG_v5.py` affecting the woven-fabric shear modulus. The remaining findings are minor or latent.

Severity levels: **CRITICAL** (non-conservative in a safety context) | **MAJOR** (theoretically incorrect or uncertified) | **MINOR** (latent bug or hardcoded assumption)

---

## ERRORS in `flutterEstimate_v3.py`

---

### ERROR F-1 — `ackeret_factor`: Wrong supersonic aerodynamic centre value
**Severity: CRITICAL**
**Location:** `flutterEstimate_v3.py:96–112`, specifically `e_sup = 0.10`

**Description:**
The Ackeret correction factor is derived from the ratio of subsonic-to-supersonic aerodynamic torsional coupling:

```python
e_sub, e_sup = 0.25, 0.10
cl_sub  = 2.0 * math.pi
cl_sup  = 4.0 / math.sqrt(M**2 - 1.0)
return math.sqrt((cl_sub * e_sub) / (cl_sup * e_sup))
```

The docstring labels `e_sub = 0.25` as "(quarter-chord AC)" and `e_sup = 0.10` as "(Ackeret AC)". Both usages are incorrect:

1. **Ackeret AC is at 0.50c, not 0.10c.** For a symmetric thin airfoil in inviscid supersonic flow, Ackeret thin-airfoil theory places the aerodynamic centre (pressure centre for a flat-plate lift distribution) at mid-chord (0.50c). The value 0.10c has no basis in standard supersonic thin-airfoil theory.

2. **AC positions are used as moment arms.** For a flutter analysis the relevant quantity is the *moment arm* from the aerodynamic centre to the fin's elastic axis, not the AC position from the leading edge. These are different quantities. The code is inconsistent: `e_sub` appears to be an AC position while `e_sup` may have been intended as a moment arm — yet the formula treats both identically.

**Consequence:**
At design point M = 1.942:
- Code produces: `f_sup = sqrt((2π × 0.25) / ((4/1.665) × 0.10)) ≈ 2.55`
- Flutter speed is inflated to `Vf_supersonic ≈ 2.55 × Vf_subsonic` — more than doubling the predicted flutter boundary.

A corrected formulation using moment arms about the elastic axis (assumed at 0.40c, consistent with typical fins) and the correct Ackeret AC (0.50c) gives:

```
e_sub_arm = 0.40 - 0.25 = 0.15c   (sub AC at 0.25c, EA at 0.40c)
e_sup_arm = 0.50 - 0.40 = 0.10c   (super AC at 0.50c, EA at 0.40c)
f_sup = sqrt((2π × 0.15) / ((4/1.665) × 0.10)) ≈ 1.98
```

Or for EA at 0.35c: `f_sup ≈ 1.32`.
The corrected factor is **22–48% lower** than the code currently computes.

**How to verify:**
- NACA TM 317 (Ackeret 1925): confirms CL_α = 4/√(M²−1) and AC at 0.50c for symmetric airfoil.
- Bisplinghoff, Ashley & Halfman (1955) §5.5: flutter moment arm derivation.
- Check: set `M = 1.0` — `ackeret_factor` should return 1.0 ✓; at M → ∞ the factor should grow as `M^(1/2)`, which this formula does for the current values — but at the wrong magnitude.

**Fix required:**
1. Define fin elastic axis position `x_ea` (fraction of chord). For a composite fin, this should come from the FEA model; a conservative default is 0.35–0.40c.
2. Compute `e_sub = x_ea - 0.25` (moment arm, subsonic).
3. Compute `e_sup = 0.50 - x_ea` (moment arm, supersonic Ackeret).
4. Ensure both `e_sub` and `e_sup` are positive (if EA is between 0.25c and 0.50c, both are positive and the formula is physically meaningful).

---

### ERROR F-2 — `ackeret_factor` applied at rocket Mach, not flutter Mach
**Severity: MAJOR**
**Location:** `flutterEstimate_v3.py:144`

```python
f_sup = ackeret_factor(M_rocket) if do_super else 1.0
```

**Description:**
The NACA 4197 subsonic formula provides a flutter speed `Vf1` and a corresponding flutter Mach `Mf1`. The ackeret correction should account for the aerodynamic properties *at the flutter condition*, i.e., should use `Mf1` (or more correctly, be solved self-consistently as a supersonic flutter problem). Using `M_rocket` (the flight Mach, which is below `Mf1` by the FSF margin) mixes the aerodynamic model evaluated at flight conditions with a structural model that defines the flutter boundary.

**Consequence:**
The correction factor is evaluated at the wrong Mach number. Since `M_rocket < Mf1` (we require FSF > 1.5), the factor computed at the lower `M_rocket` may differ from the factor at `Mf1`. For Mach numbers in the range 1.5–3, the Ackeret CL_α varies significantly. This introduces a secondary inconsistency on top of ERROR F-1.

**How to verify:**
Implement a simple iterative solver: assume a flutter Mach `Mf`, compute `ackeret_factor(Mf)`, compute `Vf = Vf_sub * f_sup * f_sw`, compute `Mf_new = Vf/a`, iterate until convergence. Compare against the non-iterative result.

**Fix required:**
Either (a) solve self-consistently for `Mf` by iteration (preferred), or (b) document explicitly that `M_rocket` is used as a conservative lower-bound for the correction factor and verify that this is indeed conservative.

---

### ERROR F-3 — `sweep_factor`: Exponent not validated against reference
**Severity: MAJOR**
**Location:** `flutterEstimate_v3.py:114–123`

```python
return 1.0 / math.sqrt(math.cos(math.radians(sweep_deg)))
```

**Description:**
The formula `V_swept = V_unswept / cos^(−1/2)(Λ)` is attributed to Bisplinghoff/Ashley/Halfman (1955) §5.5. Standard sweep correction formulations in the flutter literature use different exponents:

| Source | Formula | Exponent on cos(Λ) |
|--------|----------|--------------------|
| Bisplinghoff/Ashley/Halfman (1955) §5.5 | q_flutter ∝ cos²(Λ) | Vf ∝ cos(Λ) → factor = 1/cos(Λ) |
| Theodorsen/Garrick (NACA Rept. 685) | AR_eff ∝ cos(Λ) | varies |
| Code (current) | Vf ∝ cos^(−1/2)(Λ) | factor = cos^(−1/2)(Λ) |

At Λ = 57.4°:
- `cos(57.4°) = 0.538`
- Code factor = 1/√0.538 = **1.364** (adds 36% to flutter speed)
- BAH §5.5 cosine form = 1/cos(57.4°) = **1.858** (adds 86%)

The code produces a factor significantly lower than what the BAH formula would give, making the sweep correction **less conservative**. The source equation and derivation path must be identified and cited.

**How to verify:**
- Read Bisplinghoff/Ashley/Halfman (1955) §5.5 directly; reproduce the derivation step that leads to the cos^n(Λ) factor.
- Cross-check with NACA TN 4197 Table 2 and/or AIAA S-080 for swept fin corrections.
- Run with `--sweep-off` and compare; the sweep factor currently gives the most optimistic (lowest) penalty at the critical Λ = 57.4°.

**Fix required:**
Replace with the exponent derived from the cited reference. Add an inline equation comment with the form `# Eq. X.Y — [Reference, year]` as required by CLAUDE.md coding standards.

---

### ERROR F-4 — `flutter_naca4197`: K-factor formula is unvalidated linearisation
**Severity: MAJOR**
**Location:** `flutterEstimate_v3.py:87–88`

```python
K  = 0.65 * (1.0 + 0.10*(AR - 1.0))
K  = max(K, 0.40)
```

**Description:**
NACA TN 4197 Table 2 provides K as a function of AR and taper ratio at specific discrete points. The code implements a single-variable linear formula derived from a single reference point (K = 0.65 at AR = 1). No validation of this formula against the actual NACA 4197 table values is documented. Key concerns:

1. **No taper ratio dependency.** Table 2 of NACA TN 4197 is a function of both AR and λ = ct/cr. The code ignores λ.
2. **Linear extrapolation to low AR.** At the design AR = 0.711, K = 0.63. The accuracy of the linearisation in this low-AR regime is unconfirmed.
3. **No upper bound.** K grows unboundedly for large AR; only a lower floor (`max(K, 0.40)`) is applied.

**How to verify:**
- Obtain NACA TN 4197 Table 2 (or Figure 6 if Table 2 is not available).
- Read K values at AR = 0.5, 1.0, 2.0 for λ = 0.5 (the design taper ratio).
- Compare against the linear formula output.
- Implement a 2D bilinear interpolation over (AR, λ).

**Fix required:**
Implement a 2D interpolation table from NACA TN 4197 Table 2 with AR and λ as inputs, or add a validation comment confirming that the linear formula matches the table within ±X% over the design range.

---

### ERROR F-5 — ISA function: No altitude bound check
**Severity: MINOR (latent)**
**Location:** `flutterEstimate_v3.py:46–53`

```python
def isa(h):
    """ICAO ISA troposphere (0-11 km)."""
    T = T0 - L * h
    ...
```

For h > 11 000 m, `T` becomes negative (h = 11 000 gives T = 216.65 K, h = 12 000 gives T = 210.15 K — already below the tropopause lapse rate). At h = 44 308 m, T = 0 K (division by zero). The altitude table currently only sweeps to 4 000 m, so this is not triggered in practice, but is a latent failure mode.

**Fix required:**
Add a guard: `assert 0 <= h <= 11000, f"ISA troposphere only valid 0–11 km, got {h} m"` or implement the stratosphere (constant T = 216.65 K for 11–20 km).

---

### ERROR F-6 — Flutter formula uses D66 only; composite aeroelastic effects ignored
**Severity: MINOR (known limitation, must be documented)**
**Location:** `flutterEstimate_v3.py:139`, `compute_flutter`

```python
Geff = 12.0 * D66 / t**3
```

**Description:**
The NACA 4197 formula was derived for homogeneous isotropic metallic fins. Using G_eff = 12·D66/t³ as a surrogate for the shear modulus is valid for an isotropic or quasi-isotropic laminate where D66 alone governs torsional stiffness. For a tailored laminate (β = 20°, D16 ≠ 0), the torsional response involves D11, D22, D16, D26 in addition to D66. The bend-twist coupling term D16 (the tailoring effect from Weisshaar 1981) is explicitly **ignored** in the flutter speed calculation.

This is a conservative omission for a washout laminate (D16 < 0 at β = 20°), since washout reduces the flutter-driving moment and increases the true flutter speed. However, the omission should be stated explicitly, and a correction factor or note should be added.

**Fix required:**
Add a docstring note: "D16 bend-twist coupling is not included. For β = 20° (D16 < 0, washout), the true flutter speed is higher than computed — conservative. For β > 45° (D16 > 0, washin), this approach is non-conservative."

---

## ERRORS in `inplaneG_v5.py`

---

### ERROR I-1 — Woven fabric G12 not degraded by crimp knockdown
**Severity: MAJOR**
**Location:** `inplaneG_v5.py:99–103`

```python
E1_GA = 0.5 * (E1 + E2) * kc   # cross-ply average + crimp
return {
    ...
    'GA': (E1_GA, E1_GA, G12, 0.05, t_GA)  # G12 is NOT degraded
}
```

**Description:**
The crimp knockdown factor `kc = 0.92` (from Naik & Shembekar 1992) is applied to the extensional moduli E1 = E2 but **not** to G12 of the woven fabric. Naik & Shembekar (1992) document that fibre crimp reduces both the extensional *and* shear moduli of woven composites. For a plain weave, the G12 reduction due to crimp is typically of the same order as the E1/E2 reduction.

**Consequence:**
G12 of the GA90R woven ply is overestimated by ~8% (1/kc – 1 = 8.7%). This flows into A66 and D66 of the ABD matrix. For the 'ar1' layup where woven plies contribute ~37% of the ply count, the overestimate in D66 is approximately 3–5%.

**How to verify:**
- Naik & Shembekar (1992), J. Compos. Mater. 26(15):2196: Section on shear modulus crimp correction.
- Alternatively, measure G12 experimentally for the GA90R fabric and compare to `halpin_tsai(Gf12, Gm, Vf, xi=1.0)`.

**Fix required:**
Apply `kc` to G12 of the woven fabric:
```python
'GA': (E1_GA, E1_GA, G12 * kc, 0.05, t_GA)
```

---

### ERROR I-2 — Woven fabric nu12 = 0.05 is hardcoded without derivation
**Severity: MINOR**
**Location:** `inplaneG_v5.py:103`

```python
'GA': (E1_GA, E1_GA, G12, 0.05, t_GA)
```

**Description:**
The in-plane Poisson's ratio for the GA90R woven is set to 0.05 as a hardcoded constant. No micromechanical derivation or reference is cited. For a balanced 0/90 woven with E1 = E2, the effective ν12 can be estimated from CLT as:

```
ν12_woven ≈ (Q12_0 + Q12_90) / (Q11_0 + Q11_90)
          = 2*ν12_UD*(E2/(1−ν12²)) / (E1/(1−ν12²) + E2/(1−ν12²))
          = 2*ν12_UD*E2 / (E1 + E2)
```

Using UD values: `2 * 0.275 * 7.215 / (116.75 + 7.215) ≈ 0.032`

The actual value is lower than 0.05. The difference is small and has negligible effect on D66 (which is dominated by G12), but it affects D12 and D11 calculations.

**Fix required:**
Derive `nu12_GA` analytically or from reference. Replace the hardcoded `0.05` with a computed value and add a comment citing the source formula.

---

### ERROR I-3 — B matrix symmetry threshold may be too tight for tailored case
**Severity: MINOR**
**Location:** `inplaneG_v5.py:348–353`

```python
b_ref = abs(r['A66']) * r['t_total']
b_rel = b_max / b_ref if b_ref > 0 else b_max
ok    = b_rel < 1e-5
```

**Description:**
The symmetry check uses a relative threshold of 1e-5. For the tailored β = 20° case, the `make_symmetric` function constructs the bottom half by negating base angles (`-ang + beta`). Floating-point accumulation in the ABD integration can yield B values on the order of 1e-10 to 1e-12 N·m, which is well below the threshold — so the check passes. However, the check is comparing |B_max| against A66 × t rather than against D66 × 1/t (which has the same dimensions as B). Using `abs(r['D66'])` as the reference is more dimensionally appropriate:

```
B [N]     vs     D [N·m] / t [m]  = [N]
B [N]     vs     A [N/m] * t [m]  = [N]   ← current (correct dimensions but A66*t may be large)
```

The current check is dimensionally valid but the reference magnitude A66·t is larger than D66/t for thick laminates, making the relative threshold appear tighter than it is.

**Fix required:**
Use `b_ref = abs(r['D66']) / r['t_total']` as the denominator, or document why A66·t was chosen.

---

## Validation Checklist Against Theory

| Check | `flutterEstimate_v3.py` | `inplaneG_v5.py` | Status |
|-------|------------------------|-----------------|--------|
| ISA troposphere formula | `p = p0*(T/T0)^(g/RL)`, g/RL = 5.256 | — | ✓ Correct |
| AR formula `2b/(cr+ct)` | ✓ Consistent with NACA 4197 panel definition | — | ✓ |
| MAC formula `(2/3)cr(1+λ+λ²)/(1+λ)` | ✓ Standard trapezoidal formula | — | ✓ |
| Geff = 12·D66/t³ (flutter input) | ✓ Correct dimensional form | ✓ Same formula in ABD output | ✓ |
| CLPT Q_reduced | — | ✓ Standard Q11/Q22/Q12/Q66 | ✓ |
| Qbar rotation formulas (all 6 terms) | — | ✓ Verified vs Jones (1999) §2.9 | ✓ |
| ABD integration: A, B, D integrals | — | ✓ Standard dz, dz²/2, dz³/3 | ✓ |
| Halpin-Tsai ξ values (E2: ξ=2, G12: ξ=1) | — | ✓ Standard circular-fibre values | ✓ |
| Ackeret CL_α = 4/√(M²−1) | ✓ Correct | — | ✓ |
| Ackeret AC at 0.50c | **✗ Code uses 0.10c** (ERROR F-1) | — | **FAIL** |
| Sweep factor exponent | **? Unverified cos^(−1/2)** | — | **UNVERIFIED** |
| K factor from NACA 4197 Table 2 | **? Linear formula, no table comparison** | — | **UNVERIFIED** |
| Crimp knockdown on G12 | — | **✗ kc not applied to G12** (ERROR I-1) | **FAIL** |
| B = 0 symmetry for symmetric layup | — | ✓ Passes for all tested cases | ✓ |
| D11 = D22 for ±(45+β) NCF + woven stack | — | ✓ Verified analytically | ✓ |

---

## Impact on Flutter Safety Factor (Estimated)

Using design point: h = 1462 m, M_rocket = 1.942, D66 = 205.96 N·m, t = 5.356 mm, Λ = 57.4°.

| Correction | Effect on FSF | Direction |
|------------|---------------|-----------|
| F-1: Ackeret factor corrected (0.10→0.15c moment arm) | FSF × (1.98/2.55) = **−22%** | Non-conservative |
| F-1 (lower bound, EA at 0.35c): factor 1.32 | FSF × (1.32/2.55) = **−48%** | Non-conservative |
| F-3: Sweep exponent corrected (BAH cosine form) | FSF × (1.36/1.86) = **−27%** | Non-conservative |
| I-1: G12 degraded by kc, D66 ~3–5% lower | FSF × ~0.97–0.98 | Mildly non-conservative |

**Combined (F-1 mid + F-3):** FSF may be **30–60% lower** than the code currently reports. If the code currently shows FSF = 2.5–3.0, the corrected FSF could be in the range 1.0–2.1 — potentially below the 1.50 safety requirement.

---

## Recommended Actions (Priority Order)

1. **[CRITICAL] Fix ackeret_factor (ERROR F-1):** Define elastic axis position (obtain from FEA or use 0.40c as initial estimate), replace `e_sub = 0.25` and `e_sup = 0.10` with moment-arm formulation, and re-verify FSF against the 1.50 requirement.

2. **[CRITICAL] Re-verify FSF with corrected model:** After fixing F-1, rerun the flutter boundary analysis. If FSF < 1.50, design changes (more ±45 plies, reduced span, or increased thickness) are required before committing to the fin design.

3. **[MAJOR] Validate K factor (ERROR F-4):** Read NACA TN 4197 Table 2 and implement a 2D (AR, λ) interpolation or confirm the linear formula is within ±10% over 0.5 ≤ AR ≤ 2.0, 0.3 ≤ λ ≤ 0.7.

4. **[MAJOR] Validate sweep exponent (ERROR F-3):** Derive from Bisplinghoff/Ashley/Halfman §5.5 or replace with conservative `1/cos(Λ)` until verified.

5. **[MAJOR] Fix woven G12 crimp knockdown (ERROR I-1):** Apply `kc` to `G12` of the GA fabric in `material_props()`.

6. **[MAJOR] Address iterative Mach consistency (ERROR F-2):** Implement or document the M_rocket vs M_flutter issue.

7. **[MINOR] Add ISA altitude guard (ERROR F-5):** Add assertion `0 <= h <= 11000`.

8. **[MINOR] Derive woven nu12 (ERROR I-2):** Replace hardcoded `0.05` with calculated value.

---

*End of report.*
