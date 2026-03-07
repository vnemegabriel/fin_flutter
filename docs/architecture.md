# Architecture Guide

Developer reference: module dependencies, data flow, extension points, design decisions.

---

## 1. Module Dependency Graph

The architecture enforces a strict one-way dependency hierarchy. No module imports from a
layer above it.

```
┌─────────────────────────────────────────────────────────┐
│  ui/                                                    │
│  (screens, widgets, router, theme)                      │
└──────────────────────┬──────────────────────────────────┘
                       │ imports
┌──────────────────────▼──────────────────────────────────┐
│  services/                                              │
│  (ComputeService, PersistenceService, FileService)      │
└──────────────────────┬──────────────────────────────────┘
                       │ imports
┌──────────────────────▼──────────────────────────────────┐
│  modules/                                               │
│  flutter_analysis/ │ optimization/ │ openrocket/        │
└──────────────────────┬──────────────────────────────────┘
                       │ imports
┌──────────────────────▼──────────────────────────────────┐
│  core/aeroelastic/                                      │
│  (AeroelasticCoupler, StabilityMargin)                  │
└───────┬──────────────┬────────────────────────────────--┘
        │              │
┌───────▼──────┐  ┌────▼──────────────────────────────────┐
│ core/fea/    │  │ core/cfd/                              │
│ (mesh, elems,│  │ (VortexLattice, BiotSavart,            │
│ assembly,    │  │  PrandtlGlauert)                       │
│ eigensolvers)│  └───────────────────────────────────────┘
└───────┬──────┘
        │
┌───────▼──────────────────────────────────────────────────┐
│  core/materials/                                         │
│  (OrthotropicPly, CLTCalculator, Laminate,               │
│   MaterialDatabase, failure_criteria/)                   │
└───────────────────────────┬──────────────────────────────┘
                            │
┌───────────────────────────▼──────────────────────────────┐
│  core/models/                                            │
│  (FinGeometry, FlightCondition)                          │
└───────────────────────────┬──────────────────────────────┘
                            │
┌───────────────────────────▼──────────────────────────────┐
│  core/math/                                              │
│  (Matrix, ComplexMatrix, ComplexNumber, SparseMatrix)    │
└──────────────────────────────────────────────────────────┘
```

---

## 2. Analysis Pipeline Data-Flow

| Step | Input | Output | Class / Method |
|------|-------|--------|----------------|
| 1. Mesh generation | `FinGeometry`, `FEAConfig` | `Mesh` | `MeshGenerator.generateFinMesh()` |
| 2. Boundary conditions | `Mesh` | `BoundaryConditions` | `BoundaryConditions.clampedRoot()` |
| 3. Matrix assembly | `Mesh`, `LaminateABD`, `BoundaryConditions` | `Matrix K, M` (global) | `GlobalAssembler.assembleKirchhoff()` or `.assembleMindlin()` |
| 4. Eigensolve | `Matrix K, M`, nModes | `EigenResult` (eigenvalues, eigenvectors) | `EigenvalueSolver.solveGeneralized()` |
| 5. FEA result packaging | `EigenResult`, `Mesh` | `FEAResult` | `FEAEngine.analyze()` |
| 6. VLM panel build | `FinGeometry`, `CFDConfig` | `List<VLMPanel>` | `VortexLattice._buildPanels()` |
| 7. AIC assembly | `List<VLMPanel>`, α | `Matrix AIC` (incompressible) | `VortexLattice._buildAIC()` |
| 8. P-G correction | `Matrix AIC`, M | `Matrix AIC_c` (compressible) | `PrandtlGlauertCorrection.correctAIC()` |
| 9. VLM solve | `Matrix AIC_c`, `FlightCondition`, α | `CFDResult` (Γ, CL) | `VortexLattice.solve()` |
| 10. Interpolation matrix | `FEAResult`, `CFDResult` | `Matrix H` (nPanels × nDOF) | `AeroelasticCoupler._buildInterpolationMatrix()` |
| 11. Modal aero matrix | `Matrix H`, `Matrix AIC`, `Matrix Φ` | `Matrix Q_modal` (nModes × nModes) | `AeroelasticCoupler.buildModalAeroMatrix()` |
| 12. Flutter sweep | `K_modal`, `M_modal`, `Q_modal`, velocities, ρ | `FlutterResult` (V_F, V-g curves) | `FlutterSolver.solveUG()` |
| 13. Divergence solve | `Matrix K`, `Matrix A_s`, ρ | `DivergenceResult` (V_D) | `DivergenceSolver.solve()` |

---

## 3. Key Interfaces and Extension Points

### 3.1 Adding a New FEA Element

Implement `computeStiffness()` and `computeMass()` following the pattern in
`lib/core/fea/element/kirchhoff_element.dart`:

```dart
static Matrix computeStiffness(List<List<double>> nodeCoords, Matrix D)
static Matrix computeMass(List<List<double>> nodeCoords, double rhoTimesT)
static List<int> globalDofIndices(QuadElement element, int dofsPerNode)
```

Register the new element type in `GlobalAssembler.assembleGlobal()` by adding a branch
for the new `ElementType` enum value.

### 3.2 Adding a New Material

1. Add a `const OrthotropicPly` to `MaterialDatabase` in `lib/core/materials/material_database.dart`.
2. Add a corresponding entry to `assets/materials/material_database.json` (used for UI display).
3. No other changes required — the laminate and CLT code is fully generic.

### 3.3 Adding a New Optimizer

Implement the `Iterable<OptimizationIteration> optimize({...}) sync*` interface used by
`NelderMead` and `GeneticAlgorithm` in `lib/modules/optimization/`. The sync* generator
pattern allows the UI to stream progress without `StreamController` or `Isolate` wiring.

### 3.4 Adding a New Compressibility Correction

Add a static method to a new class in `lib/core/cfd/corrections/`:

```dart
static Matrix correct(Matrix aicIncompressible, FlightCondition condition)
```

Call it in `CFDEngine.analyze()` after the existing Prandtl-Glauert step.

---

## 4. Design Decisions

### 4.1 Pure-Dart Matrix Library

**Decision:** No `dart:ffi`, no native BLAS/LAPACK. All linear algebra is implemented
in `lib/core/math/matrix.dart` using `Float64List` (row-major).

**Rationale:** Avoids platform-specific build complexity (no `.so`/`.dylib` dependencies,
no CMake). Adequate performance for fin FEA meshes up to ~2000 DOF. The LU decomposition
with partial pivoting achieves numerical stability comparable to LAPACK routines for
well-conditioned systems.

**Trade-off:** Would be 10–100× slower than BLAS for large systems (>5000 DOF).
Acceptable for rocket fin analysis, where mesh DOF is typically 100–500.

### 4.2 Penalty Boundary Conditions

**Decision:** Clamped DOFs are enforced by adding a large diagonal value (10³⁰) to
K_ii and setting M_ii = 1, rather than eliminating the DOFs from the system (partitioning).

**Rationale:** Simpler assembly — no row/column bookkeeping, no index renumbering.
The 10³⁰ penalty exceeds the typical stiffness by 25+ orders of magnitude, effectively
enforcing zero displacement to 10-digit numerical accuracy.

**Trade-off:** Slightly increases the condition number of K. For the stiffness values
typical of CFRP fins (K_ii ~ 10³–10⁸ N/m), the ratio 10³⁰/K_ii ~ 10²²–10²⁷, well
beyond the 10¹⁶ double-precision range, so the constrained DOFs are effectively locked.

### 4.3 Quasi-Steady Aerodynamics

**Decision:** The modal aerodynamic force matrix Q_modal is real (not frequency-dependent),
using the quasi-steady assumption (reduced frequency k ≈ 0).

**Rationale:** Appropriate for preliminary design and for slowly oscillating structures
(f_flutter × chord / V ≪ 1). The quasi-steady approach gives conservative (lower) flutter
speed estimates, which is the safe direction for design.

**Trade-off:** Full unsteady doublet-lattice method (DLM) or Roger rational function
approximation would give more accurate flutter speeds for high-aspect-ratio fins or
low-speed flutter involving torsional modes. DLM is deferred to a future version.

### 4.4 Sync* Generators for Optimization

**Decision:** `NelderMead.optimize()` and `GeneticAlgorithm.optimize()` are `sync*`
generator functions returning `Iterable<OptimizationIteration>`.

**Rationale:** Allows the UI to consume iterations lazily (e.g., via a `StreamBuilder`
driven by `Stream.fromIterable()`) without needing a separate `Isolate` or
`StreamController` per optimizer. The calling code calls `.toList()` for headless use
or `.take(maxIter)` to limit iterations.

### 4.5 `Isolate.run()` for Full Analysis

**Decision:** `ComputeService.runFlutterAnalysis()` wraps `FlutterAnalysisModule.analyze()`
in `Isolate.run()`.

**Rationale:** The FEA eigensolve and VLM AIC assembly can take 1–10 s on fine meshes.
Running on the main isolate would freeze the Flutter UI thread. `Isolate.run()` is the
simplest one-shot isolate API (Dart ≥ 2.19) with automatic message passing.

---

## 5. Performance Notes

| Operation | Complexity | Typical Time |
|-----------|------------|--------------|
| LU decomposition (`Matrix.solve`) | O(n³) | n=105 DOF: < 10 ms |
| Inverse iteration (1 mode) | O(n³ × iterations) | n=105, 20 iter: ~50 ms |
| Full eigensolve (3 modes) | 3 × inverse iteration | n=105: ~150 ms |
| VLM AIC assembly | O(N² × panels) | 24 panels: < 1 ms |
| Flutter velocity sweep (40 pts × 3 modes) | O(40 × 3 × n²) | n=105: ~30 ms |
| **Default mesh total (6×4 = 24 elements, 105 DOF)** | — | **< 0.5 s** |
| **Fine mesh (12×8 = 96 elements, 351 DOF)** | — | **3–5 s** |
| **Very fine mesh (20×12 = 240 elements, 741 DOF)** | — | **15–30 s** |

The dominant cost is LU factorization within the eigenvalue iteration. The VLM and
flutter sweep contribute < 5% of total time for typical configurations.

**Profiling tip:** Run `flutter run --profile -d linux` and use Dart DevTools
(Timeline view) to identify per-step timing in the `ComputeService` isolate.
