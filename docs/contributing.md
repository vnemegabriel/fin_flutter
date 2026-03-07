# Contributing Guide

## Project Structure

```
lib/
├── core/           # Analysis engine — no Flutter dependencies
│   ├── math/       # Matrix, complex, sparse
│   ├── materials/  # CLT, plies, failure criteria
│   ├── fea/        # FEA: mesh, elements, assembly, eigensolvers
│   ├── cfd/        # VLM: Biot-Savart, AIC, corrections
│   └── aeroelastic/# Coupling, stability margins
├── modules/        # Higher-level analysis pipelines
│   ├── flutter_analysis/
│   ├── divergence/
│   ├── optimization/
│   └── openrocket/
├── models/         # App-level data models
├── services/       # Compute isolation, persistence, file handling
└── ui/             # Screens, widgets, themes, routing
```

The dependency rule: **upper layers may import lower layers; lower layers must not import
upper layers.** `core/` has no imports from `modules/` or `ui/`.

---

## Running Tests

```bash
flutter test                     # all tests
flutter test --coverage          # with lcov coverage output
flutter analyze                  # static analysis
dart format --output=none --set-exit-if-changed lib/ test/   # format check
```

All new code must pass `flutter analyze` with zero warnings.

---

## Code Style

- Follow [Dart style guide](https://dart.dev/guides/language/effective-dart/style).
- Prefer `final` and immutable data classes.
- Document public APIs with `///` doc comments.
- No `print()` statements in library code; use the test framework's `expect`.
- Format with `dart format lib/ test/` before committing.

---

## Adding a New Composite Material

1. **Add a constant** in `lib/core/materials/material_database.dart`:

```dart
static const OrthotropicPly myMaterial = OrthotropicPly(
  name: 'My Material',
  E1: 120e9, E2: 8e9, G12: 5e9, nu12: 0.3,
  rho: 1500, t: 125e-6,
  xt: 1500e6, xc: 1200e6,
  yt: 40e6, yc: 160e6, s12: 65e6,
);
```

2. **Add an entry** in `assets/materials/material_database.json` with the same values.

3. **Export** the constant from `material_database.dart` and add it to the `allMaterials`
   list so it appears in the UI dropdown.

4. **Add a test** in `test/core/materials/clt_calculator_test.dart` verifying that
   `A₁₁ = Q₁₁ × h` for a `[0]ₙ` laminate of the new material.

---

## Adding a New Failure Criterion

All failure criteria live in `lib/core/materials/failure_criteria/`.

1. Create a new file, e.g., `hashin_criterion.dart`.
2. Implement the static method pattern:

```dart
class HashinCriterion {
  static double failureIndex(
    double sigma1, double sigma2, double tau12,
    OrthotropicPly ply,
  ) {
    // ... return FI; FI = 1 at failure
  }

  static double reserveFactor(
    double sigma1, double sigma2, double tau12,
    OrthotropicPly ply,
  ) => 1.0 / failureIndex(sigma1, sigma2, tau12, ply);
}
```

3. Export from the `materials` barrel file if one exists.
4. Add test cases covering failure, margin, and zero-stress edge cases.

---

## Adding a New FEA Element Type

1. **Create** `lib/core/fea/element/my_element.dart` implementing:

```dart
class MyElement {
  /// Returns a (nDof × nDof) stiffness matrix for this element.
  static Matrix computeStiffness(List<List<double>> nodeCoords, Matrix D) { ... }

  /// Returns a (nDof × nDof) mass matrix for this element.
  static Matrix computeMass(List<List<double>> nodeCoords, double rhoTimesT) { ... }

  /// Maps element DOF index to global DOF index.
  static List<int> globalDofIndices(QuadElement element, int dofsPerNode) { ... }
}
```

2. **Add** an entry to the `ElementType` enum in `lib/core/fea/mesh/mesh.dart`.
3. **Register** in `lib/core/fea/assembly/global_assembler.dart` — add a branch to
   `assembleGlobal()` that calls your element's methods.
4. **Add tests** for stiffness symmetry, mass positivity, and the patch test.

---

## Adding a New Aerodynamic Correction

Corrections are applied to the AIC matrix before aeroelastic coupling.

1. Create `lib/core/cfd/corrections/my_correction.dart`:

```dart
class MyCorrection {
  /// Returns corrected AIC matrix given the incompressible AIC and flow parameters.
  static Matrix correct(Matrix aicIncompressible, FlightCondition condition) {
    // ...
  }
}
```

2. Call `MyCorrection.correct(aic, condition)` in `lib/core/cfd/cfd_engine.dart` after
   the Prandtl-Glauert step (or replace it, with appropriate flow-condition guards).

---

## Adding a New Optimizer

Optimizers are sync* generators yielding `OptimizationIteration` objects.

1. Create `lib/modules/optimization/my_optimizer.dart`:

```dart
class MyOptimizer {
  Iterable<OptimizationIteration> optimize({
    required double Function(List<double>) objective,
    required List<double> initialPoint,
  }) sync* {
    // ... yield OptimizationIteration(bestPoint: ..., bestValue: ..., iteration: ...);
  }
}
```

2. Add the algorithm name to the optimizer selection dropdown in the UI
   (`lib/ui/screens/optimization/optimization_screen.dart`).
3. Add convergence tests in `test/modules/optimization/`.

---

## Commit Conventions

Use imperative present tense: `Add Hashin criterion`, `Fix eigenvalue convergence`,
`Update ISA atmosphere constants`.

Format:
```
<type>: <short summary>

<optional longer explanation>
```

Types: `feat`, `fix`, `docs`, `test`, `refactor`, `perf`, `chore`.

---

## Pull Request Process

1. Fork the repository and create a feature branch from `master`.
2. Implement the feature; ensure all existing tests pass.
3. Add new tests covering your change.
4. Run `flutter analyze` and `dart format` — zero warnings/errors required.
5. Open a PR with a description of what changed, why, and how to test it.
6. Respond to review comments; squash commits before merge if requested.
