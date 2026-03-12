# Architecture & Design

This document describes the system architecture, module organization, data-flow, and design decisions for fin_flutter.

---

## 1. Overview

**fin_flutter** is a two-layer computational platform:

- **Layer 1: C++ Core** (`core/cpp/`) — High-performance header-only numerical library
  - Biot-Savart, VLM, CLT, FEA solvers
  - Aeroelastic coupling and flutter analysis
  - Compiled with CMake, Eigen3 dependency

- **Layer 2: Python Orchestration** (`core/python/`) — Analysis pipeline, CLI, visualization
  - Orchestrates C++ solvers
  - File I/O, .ork parsing, plotting
  - *Status: In development*

---

## 2. C++ Core Architecture

### 2.1 Directory Structure

```
core/cpp/
├── include/fin_flutter/
│   ├── math/
│   │   ├── types.hpp         # Eigen type aliases (Matrix3d, VectorXd, etc.)
│   │   └── vec3.hpp          # Lightweight 3D vector (dot, cross, norm)
│   │
│   ├── models/
│   │   ├── fin_geometry.hpp  # Planform: span, chords, sweep
│   │   └── flight_condition.hpp # ISA 1976, velocity, altitude
│   │
│   ├── materials/
│   │   ├── orthotropic_ply.hpp    # E1, E2, G12, nu12; Q matrix
│   │   └── clt_calculator.hpp     # [A][B][D] matrix integration
│   │
│   ├── cfd/
│   │   ├── biot_savart.hpp        # Finite/semi-infinite vortex segments
│   │   └── vortex_lattice.hpp     # AIC matrix, horseshoe vortices
│   │
│   ├── fea/
│   │   └── eigenvalue_solver.hpp  # Generalized eigenvalue [K]φ = λ[M]φ
│   │
│   ├── aeroelastic/
│   │   ├── aeroelastic_coupler.hpp  # Modal projection Q_modal = ΦᵀAICΦ
│   │   └── flutter_solver.hpp       # U-g method, divergence solver
│   │
│   └── optimization/
│       └── (placeholder for Nelder-Mead, GA)
│
├── tests/
│   ├── CMakeLists.txt
│   └── test_main.cpp         # 24 benchmark cases
│
└── CMakeLists.txt            # Find Eigen3, build library & tests
```

### 2.2 Header-Only Design

**All C++ code is in `.hpp` headers** (no separate `.cpp` implementations).

**Rationale:**
- Template-heavy (Eigen matrices, vectors)
- Faster iteration during development
- Single-file includes for users
- Avoids link-time issues with templates

**Compilation:**
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

### 2.3 Dependency Graph

```
test_main.cpp
    ↓
pipeline.hpp (orchestrator)
    ├─ flutter_solver.hpp (aeroelastic)
    ├─ aeroelastic_coupler.hpp
    ├─ vortex_lattice.hpp (CFD)
    ├─ eigenvalue_solver.hpp (FEA)
    ├─ clt_calculator.hpp (materials)
    ├─ flight_condition.hpp (models)
    ├─ fin_geometry.hpp (models)
    └─ biot_savart.hpp
        ├─ vec3.hpp (math)
        └─ types.hpp (Eigen aliases)
```

**Acyclic Design:** Each module depends only on lower-level modules; no circular dependencies.

---

## 3. Core Module Responsibilities

### 3.1 `math/` — Linear Algebra

| File | Purpose |
|------|---------|
| `types.hpp` | Eigen3 type aliases (Matrix3d, VectorXd, SparseMatrix) |
| `vec3.hpp` | Lightweight 3D vector with dot(), cross(), norm(), normalized() |

**Why separate from Eigen?** Provides a clean abstraction for small-vector operations (Biot-Savart).

### 3.2 `models/` — Problem Geometry

| File | Purpose |
|------|---------|
| `fin_geometry.hpp` | Fin planform: span, root_chord, tip_chord, sweep, thickness |
| `flight_condition.hpp` | Flight state: altitude, velocity, ISA 1976 properties |

**Design:** Immutable data structures with computed properties.

### 3.3 `materials/` — Composite Laminates

| File | Purpose |
|------|---------|
| `orthotropic_ply.hpp` | Single ply: E1, E2, G12, nu12, rho, thickness; Q and Q̄ matrices |
| `clt_calculator.hpp` | Laminate [A][B][D] matrix via integration through thickness |

**Preset materials:** AS4/3501-6, T300/5208, IM7/8552, E-glass/Epoxy.

### 3.4 `cfd/` — Aerodynamics

| File | Purpose |
|------|---------|
| `biot_savart.hpp` | Induced velocity from finite & semi-infinite vortex segments |
| `vortex_lattice.hpp` | VLM solver: panels, AIC matrix, circulation solver, lift |

**Key Fix:** Horseshoe vortex trailing segments now both have Γ=+1 (corrected from Γ=−1 for left trailing), fixing AIC sign.

### 3.5 `fea/` — Structural Dynamics

| File | Purpose |
|------|---------|
| `eigenvalue_solver.hpp` | Solve generalized eigenvalue [K]{φ} = λ[M]{φ} |

**Current Scope:** Eigenvalue extraction. Full FEA mesh/element assembly not yet implemented.

### 3.6 `aeroelastic/` — Coupling & Flutter

| File | Purpose |
|------|---------|
| `aeroelastic_coupler.hpp` | Project AIC onto FEA modes: Q_modal = ΦᵀAICΦ |
| `flutter_solver.hpp` | U-g method (flutter) and divergence solver |

---

## 4. Data Flow

### 4.1 Complete Analysis Pipeline

```
1. INPUT: Fin geometry, material laminate, flight condition
   ↓
2. CLT: Compute [A][B][D] matrices
   ↓
3. AERODYNAMICS (VLM): Circulation Γ per panel, lift CL
   ↓
4. FEA: Natural frequencies ω, mode shapes Φ
   ↓
5. AEROELASTIC COUPLING: Modal aerodynamic stiffness Q_modal
   ↓
6. FLUTTER (U-g): Velocity sweep, detect g=0 crossing → V_flutter
   ↓
7. DIVERGENCE: Min eigenvalue of [K]⁻¹·A_static → V_divergence
   ↓
8. OUTPUT: V-g diagram, flutter/divergence speeds
```

---

## 5. Design Decisions

### 5.1 Why Header-Only C++?

**Pros:** Fast iteration, template-friendly, single-file distribution, no linker issues with Eigen.
**Cons:** Compilation time if many includes.

**Decision:** Standard for Eigen-based libraries (e.g., Boost). Adopted here.

### 5.2 Why Separate Python Layer?

| Responsibility | Language | Rationale |
|---|---|---|
| Numerical solvers | C++ | Performance, Eigen integration |
| File I/O, plotting, CLI | Python | Ecosystem (matplotlib, click) |
| Orchestration | Python | Easy scripting |

### 5.3 Horseshoe Vortex Sign Convention

**Initial Bug:** Left trailing vortex had Γ=−1 → AIC diagonal positive (upwash) → negative CL.

**Fix:** Both trailing vortices have Γ=+1 → AIC diagonal negative (downwash) → positive CL.

**Status:** CL sign is now correct; magnitude requires further investigation.

---

## 6. Testing Strategy

**File:** `core/cpp/tests/test_main.cpp`

**Cases:** 24 benchmarks covering CLT, ISA, VLM

**Status:**
- Cases 1–3: ✓ All pass (CLT, ISA)
- Case 4: ⚠ CL sign correct, magnitude ~50× undersized

**Acceptance:** Compare against analytical solutions with 0.2%–5% tolerance.

---

## 7. Future Enhancements

- FEA mesh assembly and element stiffness computation
- VLM corrections (sweep, dihedral, compressibility)
- Nelder-Mead and genetic algorithm optimizers
- Python bindings (pybind11)
- Comprehensive Python orchestration layer

---

## 8. Terminology

| Term | Definition |
|------|-----------|
| **AIC** | Aerodynamic Influence Coefficient; velocity at control i from unit Γ on j |
| **Bound vortex** | Main circulation vortex around fin |
| **Trailing vortex** | Shed vorticity extending downstream |
| **Horseshoe vortex** | Complete element = bound + trailing vortices |
| **Control point** | Location (3/4-chord) where flow-tangency BC enforced |
| **g-damping** | Structural margin; g = q·Q_modal/K_modal − 1; flutter at g=0 |
| **Q matrix** | Reduced stiffness for single ply |
| **Q̄ matrix** | Transformed stiffness in laminate axes |
