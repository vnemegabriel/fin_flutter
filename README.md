# Fin Flutter Analyzer

An aerospace engineering Flutter application for fin flutter and divergence analysis of amateur and professional rockets.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Flutter](https://img.shields.io/badge/Flutter-≥3.10-02569B)](https://flutter.dev)
[![Dart](https://img.shields.io/badge/Dart-≥3.0-0175C2)](https://dart.dev)

---

## Features

- **Finite Element Analysis (FEA):** Kirchhoff DKQ and Mindlin MITC4 plate elements
- **Vortex Lattice Method (CFD):** Aerodynamic panel method with Prandtl-Glauert compressibility correction
- **Classical Lamination Theory (CLT):** Full [A][B][D] matrix computation for composite laminates
- **Aeroelastic Flutter Analysis:** U-g method for flutter speed prediction
- **Static Aeroelastic Divergence:** Minimum positive eigenvalue method
- **OpenRocket Import:** Parse .ork files to extract fin geometry and flight profiles
- **Design Optimization:** Nelder-Mead simplex and Genetic Algorithm optimizers
- **Material Database:** AS4/3501-6, T300/5208, IM7/8552, E-glass/Epoxy, Kevlar/Epoxy

---

## Quick Start

```bash
git clone https://github.com/vnemegabriel/fin_flutter.git
cd fin_flutter
flutter pub get
flutter run -d linux   # or -d chrome / -d macos / -d windows
```

See [docs/INSTALL.md](docs/INSTALL.md) for platform-specific setup instructions.

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/theory.md](docs/theory.md) | Complete mathematical theory: CLT, FEA, VLM, flutter, divergence |
| [docs/test_cases.md](docs/test_cases.md) | 9 benchmark cases with analytical solutions |
| [docs/user_guide.md](docs/user_guide.md) | Screen-by-screen UI walkthrough |
| [docs/architecture.md](docs/architecture.md) | Module dependency graph, data-flow, design decisions |
| [docs/api_reference.md](docs/api_reference.md) | Headless API usage without the Flutter UI |
| [docs/contributing.md](docs/contributing.md) | How to add materials, elements, optimizers |
| [docs/INSTALL.md](docs/INSTALL.md) | Linux, macOS, Windows, Web build instructions |

---

## Architecture

```
lib/
├── core/
│   ├── math/           # Matrix, complex arithmetic, sparse matrices
│   ├── materials/      # CLT, orthotropic plies, failure criteria
│   ├── fea/            # FEA engine: mesh, elements, assembly, eigensolvers
│   ├── cfd/            # VLM aerodynamics: Biot-Savart, AIC, corrections
│   └── aeroelastic/    # Aeroelastic coupling, stability margins
├── modules/
│   ├── flutter_analysis/  # Full flutter analysis pipeline
│   ├── divergence/        # Divergence solver module
│   ├── optimization/      # Nelder-Mead, GA optimizers
│   └── openrocket/        # .ork file parser
├── models/            # App-level data models (Project, AnalysisSession)
├── services/          # Persistence, compute isolation, file handling
└── ui/                # Screens, widgets, themes, routing
```

See [docs/architecture.md](docs/architecture.md) for the full dependency graph, pipeline
data-flow table, and design decision rationale.

---

## Analysis Pipeline

1. **Geometry Input** → Define fin planform (span, chords, sweep, thickness)
2. **Material Selection** → Choose composite ply material and build laminate stack
3. **CLT Computation** → Compute [A][B][D] matrices for the laminate
4. **Mesh Generation** → Structured quad mesh over trapezoidal fin planform
5. **FEA** → Assemble [K] and [M], solve generalized eigenvalue problem for mode shapes
6. **VLM** → Discretize fin into panels, assemble AIC matrix, solve for circulation
7. **Aeroelastic Coupling** → Project AIC onto modal coordinates → Q_modal matrix
8. **Flutter Solve** → U-g method: sweep velocity, detect g=0 crossing → V_flutter
9. **Divergence Solve** → Minimum positive eigenvalue of [K]^{-1}[A_s] → V_divergence
10. **Results** → V-g diagram, natural frequencies, stability margins

---

## Running Tests

```bash
flutter test                              # all tests
flutter test test/core/math/              # matrix library
flutter test test/core/materials/         # CLT tests
flutter test test/core/fea/               # FEA element and solver tests
flutter test test/core/cfd/               # VLM tests
flutter test test/modules/optimization/   # optimizer tests
flutter test test/modules/openrocket/     # .ork parser tests
flutter test test/integration/            # end-to-end pipeline test
```

---

## Mathematical Background

### Classical Lamination Theory (CLT)

For a laminate with N plies (midplane at z=0):

```
A_ij = Σ_k Q̄_ij^(k) * (z_k - z_{k-1})
B_ij = (1/2) Σ_k Q̄_ij^(k) * (z_k² - z_{k-1}²)
D_ij = (1/3) Σ_k Q̄_ij^(k) * (z_k³ - z_{k-1}³)
```

### Flutter Condition (U-g Method)

Flutter occurs when structural damping g crosses zero:

```
g_i(V) = q(V) · Q̃_ii / K̃_ii − 1
```

Flutter speed V_F = velocity where g_i = 0 for any mode i.

### Divergence Condition

Static aeroelastic instability when:

```
det([K] - q_D·[A_s]) = 0
q_D = λ_min([K]^{-1}·[A_s])
V_D = √(2·q_D/ρ)
```

Full derivations with references in [docs/theory.md](docs/theory.md).

---

## License

MIT License — see [LICENSE](LICENSE) for details.
