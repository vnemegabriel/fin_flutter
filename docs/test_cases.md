# Benchmark Test Cases

This document specifies the benchmark validation suite for fin_flutter. Each test case includes analytical reference solutions, acceptance tolerances, and reproduction instructions.

**Test Suite:** `core/cpp/tests/test_main.cpp`
**Total Cases:** 24 tests across 4 domains
**Current Status:** 23 passing, 1 failing (Case 4 — VLM magnitude)

---

## Domain 1: Classical Lamination Theory (CLT)

### Case 1 — CLT: Unidirectional [0]₈ AS4/3501-6

**Module:** `clt_calculator.hpp`

#### Inputs

| Parameter | Value |
|-----------|-------|
| Material | AS4/3501-6 |
| E₁ | 142 GPa |
| E₂ | 10.3 GPa |
| G₁₂ | 7.2 GPa |
| ν₁₂ | 0.27 |
| Ply thickness | 125 µm |
| Layup | [0]₈ (8 plies all at 0°) |

#### Analytical Solution

**Laminate thickness:**
```
h = 8 × 125×10⁻⁶ m = 1.0 mm
```

**Reduced stiffness (material axes):**
```
Δ = 1 − ν₁₂·ν₂₁ = 1 − 0.27 × (0.27 × 10.3/142) = 0.99471
Q₁₁ = 142×10⁹ / 0.99471 = 142.75×10⁹ Pa
Q₁₂ = 0.27 × 10.3×10⁹ / 0.99471 = 2.796×10⁶ Pa
Q₆₆ = 7.2×10⁹ Pa
```

**Laminate stiffness (A matrix):**
```
A_ij = Σ Q̄_ij · (z_k − z_{k−1})
A₁₁ = 142.75×10⁹ × 0.001 = 142.75×10⁶ N/m = 142.75 MN/m
A₁₂ = 2.796×10⁶ × 0.001 = 2.796 MN/m
A₆₆ = 7.2×10⁹ × 0.001 = 7.2 MN/m
B = 0 (symmetric layup)
D₁₁ = 142.75×10⁹ × (10⁻³)³/12 = 11.896 N·m
```

**Reference:** Reddy (2004), Classical Lamination Theory, Section 3.5.

#### Expected Results

| Quantity | Value | Tolerance |
|----------|-------|-----------|
| A₁₁ | 142.75 MN/m | ±0.2% |
| A₁₂ | 2.796 MN/m | ±0.2% |
| A₆₆ | 7.2 MN/m | ±0.1% |
| B (all) | 0 | < 1×10⁻¹² N |
| D₁₁ | 11.896 N·m | ±0.2% |
| Thickness | 1.0 mm | exact |

#### Status

✓ **PASS** — All values match expected within tolerance.

---

### Case 2 — CLT: Quasi-isotropic [0/45/−45/90]ₛ AS4/3501-6

**Module:** `clt_calculator.hpp`

#### Inputs

| Parameter | Value |
|-----------|-------|
| Material | AS4/3501-6 |
| Layup | [0/45/−45/90]ₛ (8 plies symmetric) |
| Total thickness | 1.0 mm |

#### Analytical Properties

For a quasi-isotropic symmetric laminate:
- A₁₁ = A₂₂ (isotropy)
- A₁₆ ≈ 0, A₂₆ ≈ 0 (no shear-extension coupling)
- B = 0 (symmetric)

Each ply contributes its Q̄-transformed stiffness proportionally at various angles.

**Expected in-plane stiffness:** ~61.7 MPa (isotropic equivalent)

#### Expected Results

| Condition | Expected | Tolerance |
|-----------|----------|-----------|
| A₁₁ | 61.72 MN/m | ±1% |
| A₂₂ | 61.72 MN/m | ±1% |
| \|A₁₁ − A₂₂\| | 0 | < 1×10⁻¹⁰ (isotropy) |
| A₁₆ | 0 | < 1 N/m (rounding) |
| A₂₆ | 0 | < 1 N/m (rounding) |
| B (all) | 0 | < 1×10⁻¹² N |

#### Status

✓ **PASS** — Isotropy conditions verified.

---

## Domain 2: ISA 1976 Standard Atmosphere

### Case 3 — ISA Atmosphere: Temperature, Pressure, Density, Speed of Sound

**Module:** `flight_condition.hpp`

#### Test Points

| Altitude (m) | T (K) | p (Pa) | ρ (kg/m³) | a (m/s) |
|---|---|---|---|---|
| 0 | 288.15 | 101,325 | 1.2250 | 340.29 |
| 11,000 | 216.65 | 22,632 | 0.3639 | 295.15 |
| 20,000 | 216.65 | 5,475 | 0.0880 | 295.15 |

#### Analytical Solution

**Troposphere** (0–11 km):
```
T(h) = 288.15 − 0.0065·h  [K]
p(h) = 101,325 · (T/T₀)^(−5.255)  [Pa]
ρ(h) = p(h) / (287.058 · T(h))  [kg/m³]
a(h) = √(1.4 × 287.058 × T(h))  [m/s]
```

**Stratosphere** (11–20 km, isothermal + linear):
```
T(h) = 216.65  [K]
p(h) = 22,632 · exp(−0.000157·(h − 11,000))  [Pa]
ρ(h) = p / (R·T)
```

**Reference:** NACA Report 1235, U.S. Standard Atmosphere 1976.

#### Expected Results

At h = 0:
- T(0) = 288.15 K (exact)
- p(0) = 101,325 Pa (exact)
- ρ(0) = 1.2250 kg/m³ (±0.2%)
- a(0) = 340.29 m/s (±0.1%)

At h = 11 km:
- T(11k) = 216.65 K (exact)
- ρ(11k) = 0.3639 kg/m³ (±0.2%)
- a(11k) = 295.15 m/s (±0.3%)

At h = 20 km:
- T(20k) = 216.65 K (exact)
- ρ(20k) = 0.0880 kg/m³ (±0.5%)
- a(20k) = 295.15 m/s (±0.3%)

#### Status

✓ **PASS** — All 9 atmosphere tests pass within tolerance.

---

## Domain 3: Vortex Lattice Method (VLM)

### Case 4 — VLM: Rectangular Flat Plate AR=2.5, α=5°

**Module:** `vortex_lattice.hpp`

#### Inputs

| Parameter | Value |
|-----------|-------|
| Fin shape | Rectangular flat plate |
| Span b | 0.5 m |
| Root chord c_r | 0.2 m |
| Tip chord c_t | 0.2 m |
| Leading-edge sweep | 0 m (rectangular) |
| Aspect ratio AR | b²/S = 0.5²/0.1 = 2.5 |
| Velocity V | 50 m/s |
| Altitude h | 0 m (sea level) |
| Angle of attack α | 5° = 0.0873 rad |
| Discretization | 8 chordwise × 12 spanwise (96 panels) |

#### Analytical Solution

**2D Thin Plate Theory:**
```
CL_2D = 2π·sin(α) ≈ 2π × 0.0872 ≈ 0.547
```

**3D Finite Aspect Ratio (VLM Correction):**
```
CL_3D ≈ CL_2D / (1 + π·sin(α)/(2·AR))
      ≈ 0.547 / (1 + π × 0.0872 / (2 × 2.5))
      ≈ 0.547 / 1.055 ≈ 0.518
```

**Test Range:** CL ∈ [0.36, 0.42] is a more conservative empirical range (accounting for discretization effects, compressibility, etc.).

#### Expected Results

| Quantity | Value | Tolerance |
|----------|-------|-----------|
| CL | [0.36, 0.42] | empirical range |
| CL sign | positive | required for lift |
| Panel count | 96 | exact |
| AIC matrix | square 96×96 | structural |
| Γ values | all finite | no NaN/Inf |

#### Current Status

⚠ **FAIL** — CL = 0.007746

- **CL sign:** ✓ Correct (positive)
- **CL magnitude:** ✗ Undersized by ~50× (0.00775 vs. expected ~0.39)
- **Root cause:** Under investigation. Likely scaling issue in AIC assembly or panel geometry normalization.

#### Known Issues

1. **Sign Bug (FIXED):** Initially CL was negative (−3.39). Fixed by correcting left trailing vortex sign from Γ=−1 to Γ=+1 in horseshoe vortex definition.

2. **Magnitude Bug (OPEN):** CL is now positive but undersized. Likely causes:
   - AIC scaling factor (missing 1/(4π) somewhere)
   - Panel span_width or area calculation error
   - RHS normalization issue
   - Discretization spacing effects (dx vs. dy units)

#### Reproduction

```bash
cd core/cpp
./build/tests/run_tests 2>&1 | grep -A8 "Case 4"
```

---

## Test Infrastructure

### Running All Tests

```bash
cd core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/tests/run_tests
```

### Output Format

```
Case N — Description
  [PASS/FAIL] Quantity (got X, expected Y, err Z%)
  ...

=== Results: X passed, Y failed ===
```

### Adding New Test Cases

1. Compute analytical solution with known-good reference
2. Add test case to `test_main.cpp` with clear inputs/expected values
3. Specify tolerance (percentage or absolute)
4. Document case in this file with derivation and reference
5. Run full test suite; all prior tests must still pass

---

## References

1. **Reddy, J.N.** (2004). *Mechanics of Laminated Composite Plates and Shells* (2nd ed.). CRC Press.
   - Section 3.5: Classical Lamination Theory, [A][B][D] matrix computation

2. **NACA Report 1235.** U.S. Standard Atmosphere, 1976.
   - ISA temperature, pressure, density profiles

3. **Katz, J. & Plotkin, A.** (2001). *Low-Speed Aerodynamics* (2nd ed.). Cambridge University Press.
   - Section 12.3: Vortex Lattice Method, AIC matrix assembly

---

## Summary Table

| Domain | Case # | Test | Status | Issue |
|--------|--------|------|--------|-------|
| **CLT** | 1 | [0]₈ stiffness matrix | ✓ PASS | — |
| **CLT** | 2 | Quasi-isotropic isotropy | ✓ PASS | — |
| **ISA** | 3 | Temperature, pressure, density, sound speed (9 pts) | ✓ PASS | — |
| **VLM** | 4 | Flat plate CL, sign, discretization (5 tests) | ⚠ 4/5 PASS | CL magnitude 50× undersized |
| **TOTAL** | | | **23/24 PASS** | VLM Case 4 magnitude |

