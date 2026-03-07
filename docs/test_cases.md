# Benchmark Validation Test Cases

This document provides nine benchmark cases with known analytical solutions for
validating each major analysis module. Each case specifies: inputs, derivation,
expected result, tolerance, and instructions to reproduce.

---

## Case 1 — CLT: Unidirectional [0]₈ AS4/3501-6

**Module:** Classical Lamination Theory
**File:** `test/core/materials/clt_calculator_test.dart`

### Inputs

| Parameter | Value |
|-----------|-------|
| Material | AS4/3501-6 |
| E₁ | 142 GPa |
| E₂ | 10.3 GPa |
| G₁₂ | 7.2 GPa |
| ν₁₂ | 0.27 |
| t_ply | 125 µm |
| Layup | [0]₈ (8 plies all at 0°) |

### Analytical Solution

Total thickness: h = 8 × 125×10⁻⁶ = 1.0×10⁻³ m

Reciprocal Poisson: ν₂₁ = 0.27 × 10.3/142 = 0.01958

Δ = 1 − ν₁₂ν₂₁ = 1 − 0.27 × 0.01958 = 0.99471

Q₁₁ = 142×10⁹ / 0.99471 = 142.75×10⁹ Pa

For a symmetric, balanced [0]₈ laminate:

```
A₁₁ = Q₁₁ · h = 142.75×10⁹ × 1×10⁻³ = 142.75×10⁶ N/m
B = 0   (symmetric about midplane)
D₁₁ = Q₁₁ · h³/12 = 142.75×10⁹ × (10⁻³)³/12 = 11.896 N·m
```

### Expected Results

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| A₁₁ | 142.75 MN/m | ±0.2% |
| A₁₂ | ≈ 2.78 MN/m | ±0.2% |
| A₆₆ | 7.2 MN/m | ±0.1% |
| B_ij | 0 (all) | ±1×10⁻³ N (rounding) |
| D₁₁ | 11.90 N·m | ±0.2% |

**Reference:** Reddy (2004), Table 3.3.

### Reproduce in App

1. Open **Materials** screen.
2. Select material: *AS4/3501-6*.
3. Choose preset: *Unidirectional [0]₈*.
4. View the [A], [B], [D] matrix panel; confirm A₁₁ ≈ 142.7 MN/m.

---

## Case 2 — CLT: Quasi-isotropic [0/45/−45/90]_s

**Module:** Classical Lamination Theory
**File:** `test/core/materials/clt_calculator_test.dart`

### Inputs

| Parameter | Value |
|-----------|-------|
| Material | AS4/3501-6 |
| Layup | [0/45/−45/90/90/−45/45/0] (8 plies) |
| h | 1.0 mm |

### Analytical Properties

A quasi-isotropic symmetric laminate must satisfy:
- A₁₁ = A₂₂ (isotropy in-plane — equal stiffness in all directions)
- A₁₆ = A₂₆ = 0 (no in-plane shear-extension coupling)
- B = 0 (symmetric layup)

Approximate in-plane modulus:

```
E_x ≈ E_y ≈ Q₁₁(0.375 + 0.125) + Q₁₂(0.25) + ...
```

For AS4/3501-6 eight-ply QI: E_x ≈ 55–60 GPa (depends on exact Q̄ contributions).

### Expected Results

| Condition | Expected | Tolerance |
|-----------|----------|-----------|
| A₁₁ = A₂₂ | equal | |A₁₁−A₂₂|/A₁₁ < 10⁻¹⁰ |
| A₁₆ | ≈ 0 | |A₁₆| < 1 N/m |
| A₂₆ | ≈ 0 | |A₂₆| < 1 N/m |
| B_ij | 0 | max|B_ij| < 1×10⁻³ N |

**Reference:** Reddy (2004), Table 3.7; Jones (1999), Example 4.3.

---

## Case 3 — ISA 1976 Atmosphere at Standard Altitudes

**Module:** ISA Atmosphere
**File:** `test/integration/full_pipeline_test.dart`

### Standard Checkpoints

| Altitude | T (K) | p (Pa) | ρ (kg/m³) | a (m/s) |
|----------|-------|--------|-----------|---------|
| 0 m | 288.15 | 101 325 | 1.2250 | 340.29 |
| 1 000 m | 281.65 | 89 875 | 1.1117 | 336.43 |
| 11 000 m | 216.65 | 22 700 | 0.3639 | 295.15 |
| 20 000 m | 216.65 | 5 475 | 0.0880 | 295.15 |

**Source:** NOAA-S/T 76-1562 (1976).

### Expected Results

Tolerance: ±0.1% for all quantities at listed altitudes.

### Reproduce in App

1. Open **Flight Conditions** screen.
2. Set altitude slider to 0 m. Confirm ρ = 1.225 kg/m³, a = 340.3 m/s.
3. Set altitude to 11 000 m. Confirm T = 216.65 K, a = 295.15 m/s.

---

## Case 4 — VLM: Rectangular Flat Plate, AR = 5, α = 5°

**Module:** Vortex Lattice Method
**File:** `test/core/cfd/vlm_test.dart`

### Inputs

| Parameter | Value |
|-----------|-------|
| Planform | Rectangular: span = 0.5 m, chord = 0.1 m |
| Sweep | 0 (zero sweep) |
| Panels | 8 chordwise × 12 spanwise |
| AoA (α) | 5° |
| Flight | Sea level, V = 50 m/s (M ≈ 0.15) |

### Analytical Reference

Prandtl lifting-line theory (Schlichting & Truckenbrodt 1979):

```
CL_α = 2π / (1 + 2/AR) = 2π / (1 + 2/5) = 2π / 1.4 = 4.488 /rad
CL   = CL_α · α = 4.488 × (5π/180) = 0.392
```

The VLM uses discrete panels and Weissinger's ¾-chord rule, so numerical results
will differ slightly from lifting-line theory.

### Expected Results

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| CL at α = 5° | 0.36 – 0.42 | ±10% of analytical |
| CL sign | positive | — |
| AIC matrix (N×N) | finite, non-singular | — |

### Notes

The 10% tolerance accounts for the VLM panel discretization and the Weissinger
correction vs. the continuous lifting-line approximation.

---

## Case 5 — Kirchhoff Element: Stiffness Matrix Symmetry and Patch Test

**Module:** Kirchhoff DKQ Element
**File:** `test/core/fea/kirchhoff_element_test.dart`

### Test 5a: Stiffness Matrix Symmetry

For any element geometry, K_e must be symmetric:

```
|K_e[i][j] − K_e[j][i]| < 1×10⁻¹⁰   ∀ i,j
```

This is guaranteed by the Galerkin formulation Kₑ = ∫BᵀDBdA with symmetric D.

### Test 5b: Patch Test (Constant Curvature Field)

Setup: single 0.1 m × 0.1 m square element, isotropic material (E = 70 GPa, ν = 0.3,
t = 1 mm). Apply a linear deflection field that corresponds to a constant curvature.

A correct implementation of the B matrix must reproduce a constant curvature exactly
from the nodal deflections — the patch test for completeness:

```
κ_xx = M₀ / D₁₁   (constant throughout element)
```

Expected: computed strain energy = analytical value ½ κ² D₁₁ A, to within 1%.

### Test 5c: Positive Eigenvalues of K_e (Free Element)

A free (unconstrained) element should have exactly 3 zero eigenvalues (rigid-body modes:
two translations, one rotation) and 9 strictly positive eigenvalues.

---

## Case 6 — Divergence: Rectangular Plate Analytical Solution

**Module:** Divergence Solver
**File:** `test/core/fea/divergence_solver_test.dart`

### Inputs (Simplified)

A 2×2 diagonal structural stiffness matrix and aerodynamic stiffness matrix provide
an analytically trivial case:

```
K  = [[k₁, 0],  = [[1000, 0],
      [0, k₂]]      [0, 2000]]  (N/m)

A_s = [[a₁, 0],  = [[2, 0],
       [0, a₂]]      [0, 1]]    (N/m per unit q)
```

### Analytical Solution

Eigenvalues of K⁻¹ A_s:

```
λ₁ = a₁/k₁ = 2/1000 = 0.002 m²/N   →  q_D = 500 N/m²
λ₂ = a₂/k₂ = 1/2000 = 0.0005 m²/N  →  q_D = 2000 N/m²
```

Minimum positive: q_D = 500 N/m²

### Expected Results

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| q_D (divergence dynamic pressure) | 500 N/m² | ±1% |
| V_D at ρ = 1.225 kg/m³ | √(2×500/1.225) = 28.57 m/s | ±1% |

---

## Case 7 — Nelder-Mead: Sphere Function Convergence

**Module:** Nelder-Mead Optimizer
**File:** `test/modules/optimization/nelder_mead_test.dart`

### Problem

```
f(x) = Σᵢ (2xᵢ − 1)²    minimum at xᵢ = 0.5 ∀i, f_min = 0
```

Starting point: [0.8, 0.2, 0.7, 0.3] (n = 4 dimensions)

### Expected Results

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| Best objective value | < 10⁻⁶ | — |
| Each x_i | 0.5 | ±0.01 |
| Iterations to converge | < 300 | — |

---

## Case 8 — Genetic Algorithm: Sphere Function

**Module:** Genetic Algorithm
**File:** `test/modules/optimization/genetic_algorithm_test.dart`

### Problem

Same sphere function as Case 7, but solved with the GA.

| GA parameter | Value |
|-------------|-------|
| Population size | 50 |
| Generations | 100 |
| SBX η_c | 20 |
| Polynomial mutation η_m | 20 |

### Expected Results

The GA is stochastic; tolerances are wider:

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| Best objective (gen 100) | < 0.01 | — |
| Improvement: gen 0 vs gen 100 | best_100 ≤ best_0 | — |

---

## Case 9 — OpenRocket Parser: Synthetic .ork File

**Module:** OpenRocket Importer
**File:** `test/modules/openrocket/ork_parser_test.dart`

### Synthetic Input

A minimal .ork file is a ZIP archive containing `rocket.ork` (XML). The test creates
real ZIP bytes in memory using the `archive` package:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<openrocket version="1.7">
  <rocket>
    <subcomponents>
      <stage>
        <subcomponents>
          <bodytube>
            <subcomponents>
              <trapezoidfinset>
                <name>TestFin</name>
                <rootchord>0.15</rootchord>
                <tipchord>0.08</tipchord>
                <span>0.20</span>
                <sweeplength>0.05</sweeplength>
                <thickness>0.003</thickness>
              </trapezoidfinset>
            </subcomponents>
          </bodytube>
        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>
</openrocket>
```

### Expected Results

| Field | Expected | Tolerance |
|-------|----------|-----------|
| span | 0.20 m | ±0.001 m |
| rootChord | 0.15 m | ±0.001 m |
| tipChord | 0.08 m | ±0.001 m |
| sweepLength | 0.05 m | ±0.001 m |

### Reproduce in App

1. Create a `.ork` file from OpenRocket matching the geometry above.
2. File → Open .ork in the Geometry screen.
3. Verify auto-populated values match the table.

---

## Reproducing All Cases via `flutter test`

```bash
# Run all validation cases
flutter test test/core/materials/clt_calculator_test.dart   # Cases 1–2
flutter test test/integration/full_pipeline_test.dart        # Case 3
flutter test test/core/cfd/vlm_test.dart                     # Case 4
flutter test test/core/fea/kirchhoff_element_test.dart       # Case 5
flutter test test/core/fea/divergence_solver_test.dart       # Case 6
flutter test test/modules/optimization/nelder_mead_test.dart # Case 7
flutter test test/modules/optimization/genetic_algorithm_test.dart # Case 8
flutter test test/modules/openrocket/ork_parser_test.dart    # Case 9
```
