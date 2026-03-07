# Mathematical Theory Reference

This document provides the complete mathematical background for all analysis modules in
Fin Flutter Analyzer. Each section includes: governing equations, unit conventions, key
assumptions, and a pointer to the implementing Dart source file.

---

## 1. ISA 1976 Standard Atmosphere

**Source:** `lib/core/models/flight_condition.dart`

The International Standard Atmosphere (ISA) divides the lower atmosphere into layers
with a piecewise-linear temperature profile.

### 1.1 Temperature Profile

**Troposphere** (0 ≤ h ≤ 11 000 m):

```
T(h) = T₀ − L·h
```

where T₀ = 288.15 K (sea-level standard temperature) and L = 0.0065 K/m (lapse rate).

**Tropopause / lower stratosphere** (11 000 m < h ≤ 20 000 m):

```
T(h) = T₁₁  (isothermal)
```

where T₁₁ = 216.65 K.

**Upper stratosphere** (20 000 m < h ≤ 32 000 m):

```
T(h) = T₂₀ + L₂(h − 20 000)
```

where L₂ = +0.001 K/m.

### 1.2 Pressure and Density

Troposphere (power-law form):

```
p(h) = p₀ · (T(h)/T₀)^(g₀/(R·L))
ρ(h) = p(h) / (R · T(h))
```

Isothermal layer (exponential form):

```
p(h) = p₁₁ · exp(−g₀(h − h₁₁)/(R·T₁₁))
```

**Constants:**

| Symbol | Value | Units |
|--------|-------|-------|
| p₀ | 101 325.0 | Pa |
| ρ₀ | 1.2250 | kg/m³ |
| g₀ | 9.80665 | m/s² |
| R | 287.058 | J/(kg·K) |
| γ | 1.4 | — |

### 1.3 Derived Quantities

Speed of sound:

```
a(h) = √(γ · R · T(h))
```

Dynamic pressure at velocity V:

```
q = ½ · ρ · V²
```

Mach number:

```
M = V / a(h)
```

---

## 2. Fin Geometry Definitions

**Source:** `lib/core/models/fin_geometry.dart`

### 2.1 Trapezoidal Planform

A fin is defined by four parameters: span `b`, root chord `c_r`, tip chord `c_t`,
and leading-edge sweep length `x_s`.

```
  Leading edge
  |\ ← sweep x_s →
  | \________________  ← tip chord c_t
  |  \
  |   \
  |    \
  |_____\____________  ← root chord c_r
  ↑ span b
```

### 2.2 Derived Properties

**Planform area:**
```
S = ½(c_r + c_t) · b
```

**Aspect ratio** (using mirror image across root, total span = 2b):
```
AR = (2b)² / S  =  4b² / S
```

**Leading-edge sweep angle:**
```
Λ = arctan(x_s / b)   [degrees = Λ × 180/π]
```

**Taper ratio:**
```
λ = c_t / c_r   (0 < λ ≤ 1)
```

**Mean aerodynamic chord (MAC):**
```
MAC = (2/3) · c_r · (1 + λ + λ²) / (1 + λ)
```

**Chord at spanwise position y** (0 ≤ y ≤ b):
```
c(y) = c_r − (c_r − c_t) · y/b
```

**Leading-edge x-position at y:**
```
x_LE(y) = x_s · y/b
```

---

## 3. Classical Lamination Theory (CLT)

**Source:** `lib/core/materials/clt_calculator.dart`, `lib/core/materials/orthotropic_ply.dart`

### 3.1 Single Ply — Reduced Stiffness Matrix [Q]

For an orthotropic ply with principal material axes (1 = fiber direction, 2 = transverse):

```
ν₂₁ = ν₁₂ · E₂ / E₁        (reciprocal relation)
Δ   = 1 − ν₁₂ · ν₂₁

Q₁₁ = E₁ / Δ
Q₂₂ = E₂ / Δ
Q₁₂ = ν₁₂ · E₂ / Δ
Q₆₆ = G₁₂
Q₁₆ = Q₂₆ = 0
```

### 3.2 Ply Angle Transformation [Q̄]

For a ply rotated by angle θ from the laminate x-axis (m = cos θ, n = sin θ):

```
Q̄₁₁ = Q₁₁m⁴ + 2(Q₁₂ + 2Q₆₆)m²n² + Q₂₂n⁴
Q̄₂₂ = Q₁₁n⁴ + 2(Q₁₂ + 2Q₆₆)m²n² + Q₂₂m⁴
Q̄₁₂ = (Q₁₁ + Q₂₂ − 4Q₆₆)m²n² + Q₁₂(m⁴ + n⁴)
Q̄₆₆ = (Q₁₁ + Q₂₂ − 2Q₁₂ − 2Q₆₆)m²n² + Q₆₆(m⁴ + n⁴)
Q̄₁₆ = (Q₁₁ − Q₁₂ − 2Q₆₆)m³n − (Q₂₂ − Q₁₂ − 2Q₆₆)mn³
Q̄₂₆ = (Q₁₁ − Q₁₂ − 2Q₆₆)mn³ − (Q₂₂ − Q₁₂ − 2Q₆₆)m³n
```

### 3.3 Laminate ABD Matrices

With the laminate midplane at z = 0, ply k occupying z ∈ [z_{k−1}, z_k]:

```
A_ij = Σ_k  Q̄_ij^(k) (z_k − z_{k−1})
B_ij = ½ Σ_k Q̄_ij^(k) (z_k² − z_{k-1}²)
D_ij = ⅓ Σ_k Q̄_ij^(k) (z_k³ − z_{k-1}³)
```

**[A]** — in-plane stiffness (N/m)
**[B]** — bending-extension coupling (N)
**[D]** — bending stiffness (N·m)

### 3.4 Engineering Constants from Compliance

Combined 6×6 compliance matrix [abd] = [ABD]⁻¹; the in-plane portion is [a] = A⁻¹:

```
E_x  = 1 / (a₁₁ · h)
E_y  = 1 / (a₂₂ · h)
G_xy = 1 / (a₆₆ · h)
ν_xy = −a₁₂ / a₁₁
```

where h = total laminate thickness.

### 3.5 Special Laminate Types

| Type | Condition | Significance |
|------|-----------|--------------|
| Symmetric | B = 0 | No bending-extension coupling |
| Balanced | A₁₆ = A₂₆ = 0 | No in-plane shear-extension coupling |
| Quasi-isotropic | A₁₁ = A₂₂, A₁₆ = A₂₆ = 0 | Isotropic in-plane response |

**Reference:** Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells*.

---

## 4. Failure Criteria

**Source:** `lib/core/materials/failure_criteria/tsai_wu.dart`

### 4.1 Tsai-Wu Interactive Criterion

Interaction tensors:

```
F₁  = 1/X_t − 1/X_c
F₂  = 1/Y_t − 1/Y_c
F₁₁ = 1/(X_t · X_c)
F₂₂ = 1/(Y_t · Y_c)
F₆₆ = 1/S₁₂²
F₁₂ = −½ √(F₁₁ · F₂₂)   (biaxial interaction approximation)
```

Failure index (FI = 1 at failure):

```
FI = F₁σ₁ + F₂σ₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂
```

Reserve factor: RF = 1/FI (RF < 1 → failure).

### 4.2 Tsai-Hill Criterion

```
FI = (σ₁/X)² − (σ₁σ₂/X²) + (σ₂/Y)² + (τ₁₂/S₁₂)²
```

where X = X_t if σ₁ > 0, else X_c (similarly for Y).

---

## 5. Finite Element Analysis — Kirchhoff Plate Theory

**Source:** `lib/core/fea/element/kirchhoff_element.dart`

### 5.1 Assumptions (Kirchhoff-Love)

1. Plane sections remain plane and perpendicular to the mid-surface.
2. Transverse shear deformation is neglected (thin plate: h/L ≪ 1).
3. Rotations θ_x = ∂w/∂y, θ_y = −∂w/∂x are derived from the deflection field w.

### 5.2 Degrees of Freedom

Each node carries three DOF: {w, θ_x, θ_y}. A 4-node element has 12 DOF.

### 5.3 Bending Strain-Displacement (Curvature) Matrix B

```
{κ} = {κ_xx, κ_yy, 2κ_xy}ᵀ

κ_xx = ∂θ_y/∂x
κ_yy = −∂θ_x/∂y
2κ_xy = ∂θ_y/∂y − ∂θ_x/∂x
```

### 5.4 Isoparametric Mapping

Bilinear shape functions in local coordinates (ξ, η ∈ [−1, 1]):

```
N_k(ξ,η) = ¼(1 + ξ_k·ξ)(1 + η_k·η)   for k = 1,2,3,4
```

Jacobian:

```
J = [∂x/∂ξ  ∂y/∂ξ]
    [∂x/∂η  ∂y/∂η]
```

Physical derivatives: {∂N/∂x, ∂N/∂y} = J⁻¹ {∂N/∂ξ, ∂N/∂η}.

### 5.5 Element Matrices

**Stiffness** (2×2 Gauss quadrature over the element area):

```
K_e = ∫∫ Bᵀ D B dA  ≈  Σ_{g} w_g Bᵀ(ξ_g,η_g) D B(ξ_g,η_g) det(J)
```

where **D** is the 3×3 plate bending stiffness from the D-matrix:

```
D = [[D₁₁, D₁₂, D₁₆],
     [D₁₂, D₂₂, D₂₆],
     [D₁₆, D₂₆, D₆₆]]
```

**Consistent mass** (translational DOF only, rotary inertia neglected):

```
M_e = ρ·t ∫∫ Nᵀ N dA
```

### 5.6 Mindlin MITC4 Extension

For thick plates (h/L ≥ 0.1), add transverse shear DOF. Each node carries 5 DOF:
{u, v, w, θ_x, θ_y}. Element matrix is 20×20.

Shear correction factor: κ = 5/6 (Reissner-Mindlin).

Shear stiffness contribution:

```
K_s = κ ∫∫ Bsᵀ G_s Bs dA
```

where Bs is the shear strain-displacement matrix and G_s is the transverse shear modulus matrix.

---

## 6. Generalized Eigenvalue Problem

**Source:** `lib/core/fea/solver/eigenvalue_solver.dart`

### 6.1 Problem Statement

After applying boundary conditions, the free-vibration problem is:

```
[K]{φ} = ω²[M]{φ}
```

The lowest nModes eigenvalues ω² give the squared natural angular frequencies.

### 6.2 Inverse Iteration with Shift

To find the eigenvalue nearest to shift σ:

```
([K] − σ[M]) v_{n+1} = [M] v_n
v_{n+1} ← v_{n+1} / ‖v_{n+1}‖_M
```

where ‖v‖_M = √(vᵀ[M]v) is the M-norm.

Each iteration requires one LU solve of the shifted system.

### 6.3 Deflation (M-Orthogonalization)

After converging mode i, remove its component from subsequent searches:

```
v ← v − Σ_{j<i} (φ_jᵀ M v) φ_j
```

This Gram-Schmidt deflation ensures mode i+1 is M-orthogonal to all previous modes.

### 6.4 Rayleigh Quotient

Current eigenvalue estimate (converges quadratically):

```
λ = (vᵀ [K] v) / (vᵀ [M] v)
```

Natural frequency:

```
ω_i = √λ_i   [rad/s]
f_i = ω_i / (2π)  [Hz]
```

### 6.5 Convergence

Iteration stops when:

```
|λ_{n+1} − λ_n| < tol · (1 + |λ_n|)
```

with default tol = 10⁻⁸ and maximum 500 iterations.

---

## 7. Vortex Lattice Method (VLM)

**Source:** `lib/core/cfd/vlm/vortex_lattice.dart`, `lib/core/cfd/vlm/horseshoe_vortex.dart`

### 7.1 Panel Discretization

The fin planform is divided into N_chord × N_span panels.

- **Bound vortex segment:** located at ¼-chord of each panel.
- **Control point (collocation point):** located at ¾-chord of each panel.
- **Two semi-infinite trailing vortices** extend downstream from the bound segment endpoints.

### 7.2 Biot-Savart Law

**Velocity induced by a finite vortex segment** from point A to B with circulation Γ:

```
V = (Γ/4π) · [(r₁ × r₂) / |r₁ × r₂|²] · [(r̂₁ · AB̂) − (r̂₂ · AB̂)]
```

where r₁ = P − A, r₂ = P − B, P = field point.

**Semi-infinite trailing vortex** in direction d̂:

```
V = (Γ/4π) · (r × d̂) / |r × d̂|² · (1 + r·d̂/|r|)
```

### 7.3 Aerodynamic Influence Coefficient (AIC) Matrix

Element A_ij = normal velocity at control point i induced by a unit-circulation horseshoe
vortex on panel j:

```
[AIC]{Γ} = {RHS}
```

No-penetration boundary condition at each control point:

```
RHS_i = −V_∞ · n̂_i = −V cos α · n̂_x_i  −  V sin α · n̂_z_i
```

### 7.4 Lift Calculation

Kutta-Joukowski theorem per spanwise strip:

```
L_j = ρ · V_∞ · Γ_j · Δy_j
```

Total lift coefficient:

```
CL = Σ_j L_j / (q · S_ref)
```

### 7.5 Prandtl-Glauert Compressibility Correction

**Source:** `lib/core/cfd/corrections/prandtl_glauert.dart`

Valid for M < 0.85 (subsonic). Incompressible AIC is scaled:

```
β = √(1 − M²)
AIC_compressible = AIC_incompressible / β
```

This is equivalent to stretching coordinates by 1/β in the streamwise direction.

**Warning:** The P-G correction breaks down near M = 1 (transonic). The UI displays a
warning for M ≥ 0.85.

---

## 8. Aeroelastic Coupling

**Source:** `lib/core/aeroelastic/aeroelastic_coupler.dart`

### 8.1 Displacement Interpolation Matrix H

H maps structural DOF {u_s} (deflections at FEA nodes) to aerodynamic control-point
displacements {u_a}:

```
{u_a} = [H] {u_s}
```

H is built by evaluating bilinear shape functions at each VLM control point location, using
the enclosing structural element.

### 8.2 Physical Aerodynamic Stiffness

The (quasi-steady) aerodynamic force vector at structural DOF:

```
{f_aero} = [H]ᵀ [AIC] [H] {u_s}  = [Q_phys] {u_s}
```

### 8.3 Modal Projection

Let Φ (n_DOF × n_modes) be the matrix of mass-normalized mode shapes from FEA.

Generalized aerodynamic force matrix:

```
[Q_modal] = Φᵀ [Q_phys] Φ   (n_modes × n_modes)
```

Modal structural stiffness and mass:

```
K̃ = Φᵀ K Φ = diag(ω₁², ω₂², …)    (diagonal, mass-normalized)
M̃ = Φᵀ M Φ = I                      (identity, mass-normalized)
```

---

## 9. Flutter Analysis — U-g Method

**Source:** `lib/core/fea/solver/flutter_solver.dart`

### 9.1 Quasi-Steady Aeroelastic System

In modal coordinates {η}:

```
M̃ η̈ + K̃ η = q · Q_modal · η
```

where q = ½ρV² is dynamic pressure.

### 9.2 Damping Indicator

Rearranging as a per-mode stability measure:

```
g_i(V) = q(V) · Q̃_ii / K̃_ii − 1
```

Physical interpretation:
- g_i < 0: aerodynamic load is less than structural stiffness → stable
- g_i = 0: flutter onset
- g_i > 0: unstable (flutter has occurred)

### 9.3 Flutter Speed Detection

Sweep velocity V from V_min to V_max in N_steps. For each mode i:

1. Compute g_i(V_k) and g_i(V_{k+1}).
2. If g_i crosses zero (sign change), linearly interpolate:
   ```
   V_F = V_k − g_i(V_k) · (V_{k+1} − V_k) / (g_i(V_{k+1}) − g_i(V_k))
   ```
3. V_flutter = min{V_F} over all modes.

### 9.4 Output

- **V-g curves:** g_i(V) for each mode (for display in V-g diagram)
- **V_flutter:** critical flutter speed [m/s]
- **Flutter frequency:** natural frequency of the fluttering mode at V_F [Hz]

---

## 10. Divergence Analysis

**Source:** `lib/core/fea/solver/divergence_solver.dart`

### 10.1 Physical Mechanism

Divergence is a static aeroelastic instability where the aerodynamic torsional moment
increases faster with twist angle than the structural restoring moment, leading to
unbounded deflection.

### 10.2 Eigenvalue Formulation

At divergence dynamic pressure q_D:

```
[K]{x} = q_D · [A_s]{x}
```

Equivalent to:

```
q_D = minimum positive eigenvalue of  [K]⁻¹ [A_s]
```

### 10.3 Inverse Power Iteration

```
v_{n+1} = [K]⁻¹ [A_s] v_n
v_{n+1} ← v_{n+1} / ‖v_{n+1}‖

q_D ≈ vᵀ [A_s] v / vᵀ [K] v   (Rayleigh quotient)
```

The method converges to the eigenvalue of smallest magnitude. Post-processing checks
that q_D > 0 (physical requirement).

### 10.4 Divergence Speed

```
V_D = √(2 q_D / ρ)   [m/s]
```

---

## 11. Optimization

**Source:** `lib/modules/optimization/nelder_mead.dart`, `lib/modules/optimization/genetic_algorithm.dart`

### 11.1 Design Variables

Normalized to [0, 1]:

```
x_norm = (x − x_min) / (x_max − x_min)
```

Typical design vector: `[span, c_root, c_tip, sweep, θ₁, t₁, θ₂, t₂, …, θ_n, t_n]`

### 11.2 Loss Function

**Source:** `lib/modules/optimization/loss_function.dart`

```
L = w₁(V_F/V_ref − 1)² + w₂(V_D/V_ref − 1)² + w₃(m/m_ref) + μ Σ_i max(0, g_i(x))²
```

| Term | Description |
|------|-------------|
| w₁ term | Flutter margin objective |
| w₂ term | Divergence margin objective |
| w₃ term | Mass penalty |
| μ term | Constraint violation penalty |

Default weights: w₁ = 1, w₂ = 1, w₃ = 0.1, μ = 100.

### 11.3 Nelder-Mead Simplex

Operates on a simplex of n+1 vertices in n-dimensional normalized space.

| Operation | Parameter | Value |
|-----------|-----------|-------|
| Reflection | α | 1.0 |
| Expansion | γ | 2.0 |
| Contraction | ρ | 0.5 |
| Shrink | σ | 0.5 |

Algorithm:
1. Sort vertices by objective: f(x₁) ≤ f(x₂) ≤ … ≤ f(x_{n+1})
2. Compute centroid x₀ = (1/n) Σ_{i=1}^{n} x_i (excluding worst)
3. Reflect: x_r = x₀ + α(x₀ − x_{n+1})
4. If f(x_r) < f(x₁): expand; else if f(x_r) < f(x_n): accept; else: contract or shrink
5. Convergence: std dev of f values < tol, default tol = 10⁻⁶

### 11.4 Genetic Algorithm

| Parameter | Value |
|-----------|-------|
| Population size | 50 (default) |
| Crossover: SBX distribution index η_c | 20 |
| Mutation: polynomial distribution index η_m | 20 |
| Selection: tournament size k | 3 |
| Elitism count | 2 |

**SBX crossover** for parent variables u₁, u₂:

```
β = {(2r)^(1/(η_c+1))         if r ≤ 0.5
    {(1/(2(1−r)))^(1/(η_c+1))  otherwise

offspring₁ = ½[(1+β)u₁ + (1−β)u₂]
offspring₂ = ½[(1−β)u₁ + (1+β)u₂]
```

**Polynomial mutation** with mutation rate 1/n:

```
δ = {(2r)^(1/(η_m+1)) − 1           if r < 0.5
    {1 − (2(1−r))^(1/(η_m+1))        otherwise

x_mutated = x + δ · (x_max − x_min)
```

---

## 12. Stability Margins

**Source:** `lib/core/aeroelastic/stability_margin.dart`

### 12.1 Flutter Margin

```
MF = V_F / V_max
```

### 12.2 Divergence Margin

```
MD = V_D / V_max
```

where V_max is the maximum expected flight velocity.

### 12.3 Safety Status Thresholds

| Status | Condition | Color |
|--------|-----------|-------|
| SAFE | margin ≥ 1.5 | Green |
| WARNING | 1.2 ≤ margin < 1.5 | Yellow |
| MARGINAL | 1.0 ≤ margin < 1.2 | Orange |
| CRITICAL | margin < 1.0 | Red |

**Regulatory reference:** MIL-SPEC-9490D requires flutter margin ≥ 1.15 (military minimum).
The SAFE threshold of 1.5 follows common amateur rocketry design practice.

---

## References

1. Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*. CRC Press.
2. Bisplinghoff, R.L., Ashley, H., Halfman, R.L. (1955). *Aeroelasticity*. Addison-Wesley.
3. Katz, J., Plotkin, A. (2001). *Low-Speed Aerodynamics*. Cambridge University Press.
4. Bathe, K.J. (1996). *Finite Element Procedures*. Prentice Hall.
5. NOAA, NASA, USAF (1976). *U.S. Standard Atmosphere, 1976*. NOAA-S/T 76-1562.
6. Deb, K., Agrawal, R.B. (1995). Simulated binary crossover for continuous search space. *Complex Systems* 9(2), 115–148.
7. Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. *Computer Journal* 7(4), 308–313.
8. Schlichting, H., Truckenbrodt, E. (1979). *Aerodynamics of the Airplane*. McGraw-Hill.
