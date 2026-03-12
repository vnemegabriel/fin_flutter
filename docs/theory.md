# Mathematical Theory Reference

This document provides the complete mathematical background for fin_flutter's computational modules. Each section includes governing equations, unit conventions, key assumptions, and implementation references (C++ headers or Python modules).

---

## 1. ISA 1976 Standard Atmosphere

**Implementation:** `core/cpp/include/fin_flutter/models/flight_condition.hpp`

The International Standard Atmosphere (ISA) divides the lower atmosphere into layers with piecewise-linear temperature profiles.

### 1.1 Temperature Profile

**Troposphere** (0 ≤ h ≤ 11,000 m):
```
T(h) = T₀ − L·h
```
where T₀ = 288.15 K (sea-level), L = 0.0065 K/m (lapse rate).

**Tropopause** (11,000 m < h ≤ 20,000 m):
```
T(h) = 216.65 K (isothermal)
```

**Stratosphere** (20,000 m < h ≤ 32,000 m):
```
T(h) = 216.65 + 0.001·(h − 20,000) K
```

### 1.2 Pressure and Density

**Troposphere** (power law):
```
p(h) = p₀ · (T(h)/T₀)^(−g₀/(R·L))
ρ(h) = p(h) / (R·T(h))
```

**Isothermal/Stratosphere** (exponential):
```
p(h) = p_ref · exp(−g₀·(h − h_ref)/(R·T_ref))
```

**Constants:**
| Symbol | Value | Units |
|--------|-------|-------|
| p₀ | 101,325 | Pa |
| ρ₀ | 1.2250 | kg/m³ |
| g₀ | 9.80665 | m/s² |
| R | 287.058 | J/(kg·K) |
| γ | 1.4 | — |

### 1.3 Derived Quantities

Speed of sound:
```
a(h) = √(γ·R·T(h))
```

Dynamic pressure:
```
q = ½·ρ·V²
```

**Reference:** NACA Report 1235 (U.S. Standard Atmosphere, 1976)

---

## 2. Classical Lamination Theory (CLT)

**Implementation:** `core/cpp/include/fin_flutter/materials/clt_calculator.hpp`

### 2.1 Ply-Level Constitutive Relations

For an orthotropic ply with fibers at angle θ to the x-axis:

**Reduced stiffness (material axes):**
```
Q₁₁ = E₁ / (1 − ν₁₂·ν₂₁)
Q₂₂ = E₂ / (1 − ν₁₂·ν₂₁)
Q₁₂ = ν₁₂·E₂ / (1 − ν₁₂·ν₂₁)
Q₆₆ = G₁₂
```

**Transformed stiffness (laminate axes, Eq. 1.3.81, Reddy 2004):**
```
Q̄ᵢⱼ(θ) = Tᵢₘ(θ) · Qₘₙ · Tₙⱼ(θ)
```

where transformation matrix T depends on cos(θ) and sin(θ).

### 2.2 Laminate Stiffness Matrix

For a laminate with N plies stacked symmetrically about z = 0:

**Extensional stiffness:**
```
A_ij = Σₖ Q̄_ij^(k) · (z_k − z_{k−1})
```

**Bending-extension coupling:**
```
B_ij = ½ · Σₖ Q̄_ij^(k) · (z_k² − z_{k−1}²)
```

**Bending stiffness:**
```
D_ij = ⅓ · Σₖ Q̄_ij^(k) · (z_k³ − z_{k−1}³)
```

where z_k is the coordinate of the k-th ply interface.

**Key Assumptions:**
- Perfect bonding between plies (no slip)
- Linear elastic material behavior
- Small strains and rotations
- Plane-stress assumption in each ply

**Reference:** Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells* (2nd ed.), Eqs. 3.5–3.7.

---

## 3. Vortex Lattice Method (VLM)

**Implementation:** `core/cpp/include/fin_flutter/cfd/vortex_lattice.hpp`

### 3.1 Horseshoe Vortex Element

Each panel carries a horseshoe vortex consisting of:

1. **Bound vortex** (at 1/4-chord): filament from A (inboard, y₀) to B (outboard, y₁) with strength Γ
2. **Right trailing vortex** (from B): extends to +∞ downstream in +x with strength Γ
3. **Left trailing vortex** (from A): represents the return leg, effectively Γ at A extending downstream

**Control point location:** 3/4-chord (Pistolesi theorem for flat-plate approximation)

### 3.2 Biot-Savart Law for Induced Velocity

**Finite segment** (bound vortex from A to B):
```
V = (Γ/4π) · (r₁ × r₂) / |r₁ × r₂|² · (cos_term / |AB|)
```

where cos_term = AB·(r̂₁ − r̂₂), r₁ = P − A, r₂ = P − B.

**Semi-infinite segment** (trailing vortex from A in direction d̂):
```
V = (Γ/4π) · (r × d̂) / |r × d̂|² · (1 + r·d̂/|r|)
```

**Reference:** Katz, J. & Plotkin, A. (2001). *Low-Speed Aerodynamics* (2nd ed.), Eqs. 10.11, 10.15.

### 3.3 Aerodynamic Influence Coefficient (AIC) Matrix

AIC[i,j] = normal velocity at control point i induced by unit circulation on panel j:
```
AIC_ij = n̂ᵢ · (v_bound + v_trail_A + v_trail_B)
```

**Flow-tangency boundary condition:**
```
(V_freestream + Σⱼ AIC_ij · Γⱼ) · n̂ᵢ = 0
```

Rearranged:
```
AIC · Γ = −V_freestream · n̂
```

### 3.4 Lift Coefficient

**Kutta-Joukowski theorem (per panel):**
```
L_i = ρ·V_∞·Γᵢ·Δyᵢ
```

**Total lift coefficient:**
```
CL = (Σᵢ L_i) / (q·S)
```

where q = ½·ρ·V∞², S = planform area.

### 3.5 Current Status & Known Issues

**Implemented:** ✓ Horseshoe vortex discretization, AIC assembly, Γ solver, lift computation

**Known Issue:** Test Case 4 (flat plate, AR=2.5, α=5°) produces CL = 0.0077 (positive sign correct, but undersized by ~50×). Under investigation.

---

## 4. Finite Element Analysis (FEA)

**Implementation:** `core/cpp/include/fin_flutter/fea/eigenvalue_solver.hpp`

### 4.1 Generalized Eigenvalue Problem

Modal analysis of a fin structure:
```
[K]{φᵢ} = λᵢ[M]{φᵢ}
```

where:
- [K] = stiffness matrix (from CLT and laminate geometry)
- [M] = consistent mass matrix
- λᵢ = (ωᵢ)² (eigenvalue = natural frequency squared)
- {φᵢ} = mass-normalized eigenvector

**Solver:** `Eigen::GeneralizedSelfAdjointEigenSolver` (QZ decomposition)

**Output:** Eigenvalues sorted ascending, eigenvectors M-normalized.

---

## 5. Aeroelastic Coupling

**Implementation:** `core/cpp/include/fin_flutter/aeroelastic/aeroelastic_coupler.hpp`

### 5.1 Modal Aerodynamic Matrix

Project the AIC matrix onto the FEA modal basis:
```
Q_modal_ij = Φᵀ · Hᵀ · AIC · H · Φ
```

where:
- Φ = [φ₁ | φ₂ | ... | φₙ] (modal matrix, N × n_modes)
- H = interpolation matrix mapping panel collocation points to mesh nodes
- AIC = full aerodynamic influence coefficient matrix (N_panels × N_panels)

**H assembly:** Nearest-node spanwise interpolation.

---

## 6. Flutter Analysis (U-g Method)

**Implementation:** `core/cpp/include/fin_flutter/aeroelastic/flutter_solver.hpp`

### 6.1 Structural Damping Parameter

For each mode i:
```
g_i(V) = q(V) · Q̃_modal_ii / K̃_modal_ii − 1
```

where:
- q(V) = ½·ρ(h)·V² (dynamic pressure)
- Q̃_modal_ii = modal aerodynamic stiffness
- K̃_modal_ii = modal structural stiffness

**Flutter condition:** g_i = 0

**Method:** Sweep velocity V from 0 to V_max, detect zero-crossing of g.

### 6.2 Output

VGPoint structure:
- Velocity V [m/s]
- Damping g [—] for each mode
- Flutter margin (V_flutter − V) / V

**Reference:** Theodorsen, T. & Garrick, I.E. (1942). Mechanism of flutter — A theoretical and experimental investigation. NACA Report 685.

---

## 7. Static Divergence Analysis

**Implementation:** `core/cpp/include/fin_flutter/aeroelastic/flutter_solver.hpp` (DivergenceSolver)

### 7.1 Divergence Eigenvalue Problem

Static aeroelastic instability occurs when:
```
det([K] − q_D·[A_static]) = 0
```

where [A_static] = Hᵀ · AIC · H (static aerodynamic stiffness matrix in structural coordinates).

**Solution:** Minimum positive eigenvalue λ_min([K]⁻¹·[A_static])

**Divergence dynamic pressure:**
```
q_D = 1 / λ_min
V_D = √(2·q_D / ρ)
```

---

## 8. Material Database

**Storage:** `assets/materials.json` (planned)

### 8.1 Implemented Materials

| Material | E₁ (GPa) | E₂ (GPa) | G₁₂ (GPa) | ν₁₂ | ρ (kg/m³) | Source |
|----------|----------|----------|-----------|-----|-----------|--------|
| AS4/3501-6 | 142 | 10.3 | 7.2 | 0.27 | 1600 | MTM35 datasheet |
| T300/5208 | 138 | 8.96 | 7.1 | 0.30 | 1600 | Hexcel |
| IM7/8552 | 171 | 10.3 | 7.4 | 0.27 | 1600 | Hexcel |
| E-glass/Epoxy | 45.6 | 12.0 | 4.5 | 0.27 | 1950 | Generic |

---

## 9. Units & Conventions

**SI Units Throughout:**
- Length: meters [m]
- Force: Newtons [N]
- Pressure/Stress: Pascals [Pa]
- Mass: kilograms [kg]
- Temperature: Kelvin [K]
- Angle: radians [rad]

**Coordinate System:**
- **x**: chordwise direction (leading edge → trailing edge)
- **y**: spanwise direction (root → tip)
- **z**: normal to the plate (pointing "up" for lift)

**Laminate Stacking Convention:**
- Plies indexed from root (z = 0) to tip
- Symmetric laminates: [0/45/−45/90]_s means 8 plies total

---

## 10. References

1. Bisplinghoff, R.L., Ashley, H. & Halfman, R.L. (1955). *Aeroelasticity*. Addison-Wesley.
2. Katz, J. & Plotkin, A. (2001). *Low-Speed Aerodynamics* (2nd ed.). Cambridge University Press. [Sections 12.1–12.3]
3. Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells* (2nd ed.). CRC Press. [Sections 3.1–3.2]
4. NASA Technical Report Server (NTRS). ISA 1976 Standard Atmosphere. NACA Report 1235.
5. Theodorsen, T. & Garrick, I.E. (1942). Mechanism of flutter. NACA Report 685.

See `Recursos/` for PDFs of key papers.
