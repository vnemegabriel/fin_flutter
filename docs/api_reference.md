# API Reference — Headless Use

The analysis engine can be used programmatically without the Flutter UI.
All core classes are in `lib/core/` and `lib/modules/`.

---

## Quick Start

```dart
import 'package:fin_flutter/core/models/fin_geometry.dart';
import 'package:fin_flutter/core/models/flight_condition.dart';
import 'package:fin_flutter/core/materials/material_database.dart';
import 'package:fin_flutter/core/materials/laminate.dart';
import 'package:fin_flutter/modules/flutter_analysis/flutter_analysis_module.dart';
import 'package:fin_flutter/modules/flutter_analysis/flutter_input.dart';
import 'package:fin_flutter/core/fea/fea_engine.dart';
import 'package:fin_flutter/core/cfd/cfd_engine.dart';

void main() {
  // 1. Define geometry
  final geometry = FinGeometry(
    span: 0.15,
    rootChord: 0.12,
    tipChord: 0.06,
    sweepLength: 0.04,
    thickness: 0.003,
  );

  // 2. Build laminate
  final ply = MaterialDatabase.im78552;
  final laminate = Laminate.quasiIsotropic(ply);
  final abd = laminate.computeABD();

  // 3. Set flight condition
  final condition = FlightCondition.isa(altitude: 1000, mach: 0.5);

  // 4. Run full analysis
  final module = FlutterAnalysisModule();
  final output = module.analyze(FlutterAnalysisInput(
    geometry: geometry,
    laminate: abd,
    flightCondition: condition,
    feaConfig: const FEAConfig(
      spanwiseElements: 6,
      chordwiseElements: 4,
      elementType: ElementType.kirchhoffDKQ,
      nModes: 3,
    ),
    cfdConfig: const CFDConfig(
      chordwisePanels: 4,
      spanwisePanels: 6,
      angleOfAttackDeg: 2.0,
    ),
    velocityMin: 20.0,
    velocityMax: 500.0,
    velocitySteps: 40,
    maxFlightVelocity: 200.0,
  ));

  print('Flutter speed:   ${output.flutterResult.flutterSpeed?.toStringAsFixed(1)} m/s');
  print('Divergence speed:${output.divergenceResult?.divergenceSpeed.toStringAsFixed(1)} m/s');
  print('Flutter margin:  ${output.flutterMargin?.toStringAsFixed(3)}');
}
```

---

## Core Classes

### `FinGeometry`

`lib/core/models/fin_geometry.dart`

```dart
FinGeometry({
  required double span,          // [m] semi-span (one fin)
  required double rootChord,     // [m]
  required double tipChord,      // [m]
  required double sweepLength,   // [m] leading-edge offset at tip
  required double thickness,     // [m] total fin thickness
})
```

**Computed properties:**

| Property | Type | Description |
|----------|------|-------------|
| `planformArea` | `double` | S = ½(c_r + c_t) b  [m²] |
| `aspectRatio` | `double` | AR = 4b² / S |
| `sweepAngleDeg` | `double` | Λ = arctan(x_s/b) [°] |
| `taperRatio` | `double` | λ = c_t / c_r |
| `chordAt(double y)` | `double` | Chord at spanwise station y [m] |
| `leadingEdgeX(double y)` | `double` | LE x-position at y [m] |

---

### `FlightCondition`

`lib/core/models/flight_condition.dart`

```dart
// ISA atmosphere factory (most common usage)
FlightCondition.isa({required double altitude, required double mach})

// Direct construction
FlightCondition({
  required double density,       // ρ [kg/m³]
  required double speedOfSound,  // a [m/s]
  required double velocity,      // V [m/s]
  required double mach,          // M = V/a
})
```

**Properties:** `density`, `velocity`, `mach`, `speedOfSound`, `dynamicPressure`.

---

### `OrthotropicPly`

`lib/core/materials/orthotropic_ply.dart`

```dart
OrthotropicPly({
  required String name,
  required double E1,    // Fiber-direction modulus [Pa]
  required double E2,    // Transverse modulus [Pa]
  required double G12,   // In-plane shear modulus [Pa]
  required double nu12,  // Major Poisson ratio
  required double rho,   // Density [kg/m³]
  required double t,     // Nominal ply thickness [m]
  // Strength (optional, for failure criteria):
  double xt = 0, double xc = 0,  // Fiber tension/compression strength [Pa]
  double yt = 0, double yc = 0,  // Transverse tension/compression strength [Pa]
  double s12 = 0,                 // In-plane shear strength [Pa]
})
```

**Methods:**
- `reducedStiffness` → `Matrix` (3×3 [Q] matrix)
- `transformedStiffness(double thetaDeg)` → `Matrix` (3×3 [Q̄] at angle θ)

---

### `MaterialDatabase`

`lib/core/materials/material_database.dart`

Pre-defined constants:

```dart
MaterialDatabase.as43501     // AS4/3501-6 carbon/epoxy
MaterialDatabase.t3005208    // T300/5208 carbon/epoxy
MaterialDatabase.im78552     // IM7/8552 carbon/epoxy
MaterialDatabase.eglassEpoxy // E-glass/Epoxy
MaterialDatabase.kevlarEpoxy // Kevlar/Epoxy
MaterialDatabase.cfrpGeneric // Generic CFRP (approximate)
```

---

### `Laminate`

`lib/core/materials/laminate.dart`

```dart
// Factory constructors
Laminate.unidirectional(OrthotropicPly ply, {int nPlies = 8})
Laminate.crossPly(OrthotropicPly ply, int pairs)        // [0/90]_pairs_s
Laminate.quasiIsotropic(OrthotropicPly ply)             // [0/45/-45/90]_s

// Manual construction
final laminate = Laminate(plies: [
  LaminatePly(ply: material, angleDeg: 0),
  LaminatePly(ply: material, angleDeg: 45),
  LaminatePly(ply: material, angleDeg: -45),
  LaminatePly(ply: material, angleDeg: 90),
  LaminatePly(ply: material, angleDeg: 90),
  LaminatePly(ply: material, angleDeg: -45),
  LaminatePly(ply: material, angleDeg: 45),
  LaminatePly(ply: material, angleDeg: 0),
]);
```

**Method:**
- `computeABD()` → `LaminateABD`

---

### `LaminateABD`

Result of `Laminate.computeABD()`:

```dart
LaminateABD {
  Matrix A;          // In-plane stiffness (3×3) [N/m]
  Matrix B;          // Bending-extension coupling (3×3) [N]
  Matrix D;          // Bending stiffness (3×3) [N·m]
  Matrix combinedABD;// Full 6×6 matrix
  double thickness;  // Total laminate thickness [m]
}
```

---

### `FEAEngine`

`lib/core/fea/fea_engine.dart`

```dart
final engine = FEAEngine(
  config: FEAConfig(
    spanwiseElements: 6,
    chordwiseElements: 4,
    elementType: ElementType.kirchhoffDKQ,  // or mindlinMITC4
    nModes: 3,
  ),
);

final FEAResult result = engine.analyze(
  geometry: geometry,
  abd: abd,
);
```

**`FEAResult` fields:**

| Field | Type | Description |
|-------|------|-------------|
| `naturalFrequenciesHz` | `List<double>` | Sorted natural frequencies [Hz] |
| `modeShapes` | `List<List<double>>` | Modal displacement vectors |
| `nModes` | `int` | Number of converged modes |
| `mesh` | `Mesh` | The finite element mesh |
| `K` | `Matrix` | Global stiffness matrix |
| `M` | `Matrix` | Global mass matrix |

---

### `CFDEngine`

`lib/core/cfd/cfd_engine.dart`

```dart
final engine = CFDEngine(
  config: const CFDConfig(
    chordwisePanels: 4,
    spanwisePanels: 6,
    angleOfAttackDeg: 2.0,
  ),
);

final CFDResult result = engine.analyze(
  geometry: geometry,
  condition: condition,
);
```

**`CFDResult` fields:** `CL`, `circulationStrengths`, `controlPoints`, `aic`.

---

### `FlutterAnalysisModule`

`lib/modules/flutter_analysis/flutter_analysis_module.dart`

Orchestrates FEA → VLM → Coupling → Flutter → Divergence.

```dart
final FlutterAnalysisOutput output = FlutterAnalysisModule().analyze(input);
```

**`FlutterAnalysisOutput` fields:**

| Field | Type | Description |
|-------|------|-------------|
| `feaResult` | `FEAResult` | Structural modes |
| `cfdResult` | `CFDResult` | Aerodynamic solution |
| `flutterResult` | `FlutterResult` | V_flutter, V-g curves |
| `divergenceResult` | `DivergenceResult?` | V_divergence (null if none found) |
| `flutterMargin` | `double?` | V_F / V_max |
| `divergenceMargin` | `double?` | V_D / V_max |
| `stabilityStatus` | `SafetyStatus` | SAFE / WARNING / MARGINAL / CRITICAL |

---

### `OrkImporter`

`lib/modules/openrocket/ork_parser.dart`

```dart
import 'dart:io';

final bytes = await File('rocket.ork').readAsBytes();
final OrkImportResult result = OrkParser.parse(bytes);

final FinGeometry geometry = result.finGeometry;
// result.flightData  →  List<FlightDataPoint> (altitude, velocity, time)
```

Throws `OrkParseException` for invalid or corrupt files.

---

### `NelderMead`

`lib/modules/optimization/nelder_mead.dart`

```dart
final optimizer = NelderMead(
  maxIterations: 500,
  convergenceTolerance: 1e-6,
);

for (final iteration in optimizer.optimize(
  objective: (x) => myLossFunction(x),
  initialPoint: [0.5, 0.5, 0.5],
)) {
  print('Iter ${iteration.iteration}: ${iteration.bestValue}');
}
```

Each yielded `OptimizationIteration` has: `bestPoint`, `bestValue`, `iteration`.

---

### `GeneticAlgorithm`

`lib/modules/optimization/genetic_algorithm.dart`

```dart
final ga = GeneticAlgorithm(
  populationSize: 50,
  maxGenerations: 200,
  crossoverEta: 20,
  mutationEta: 20,
);

for (final iter in ga.optimize(
  objective: (x) => myLossFunction(x),
  dimensions: 6,
)) {
  // iter.bestValue, iter.bestPoint, iter.iteration (generation number)
}
```

---

### `ComputeService` (Isolate wrapper)

`lib/services/compute_service.dart`

For use in Flutter apps — dispatches analysis to a background isolate:

```dart
final output = await ComputeService.runFlutterAnalysis(input);
```

Returns the same `FlutterAnalysisOutput` as `FlutterAnalysisModule.analyze()`, but
runs in a separate Dart isolate to keep the UI thread responsive.

---

## Error Handling

| Exception | Thrown by | Cause |
|-----------|-----------|-------|
| `OrkParseException` | `OrkParser.parse()` | Invalid ZIP or missing XML elements |
| `EigensolverException` | `EigenvalueSolver` | Non-convergence after max iterations |
| `MatrixSingularException` | `Matrix.solve()` | Singular or near-singular matrix |
| `StateError` | Mesh generator | Zero-area element (degenerate geometry) |

All exceptions carry a descriptive message. Wrap calls in try/catch for production use.
