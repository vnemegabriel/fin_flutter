# Fin Flutter Analyzer

An aerospace engineering Flutter application for fin flutter and divergence analysis of amateur and professional rockets.

## Features

- **Finite Element Analysis (FEA)**: Kirchhoff DKQ and Mindlin MITC4 plate elements
- **Vortex Lattice Method (CFD)**: Aerodynamic panel method with Prandtl-Glauert compressibility correction
- **Classical Lamination Theory (CLT)**: Full [A][B][D] matrix computation for composite laminates
- **Aeroelastic Flutter Analysis**: U-g method for flutter speed prediction
- **Static Aeroelastic Divergence**: Minimum positive eigenvalue method
- **OpenRocket Import**: Parse .ork files to extract fin geometry and flight profiles
- **Design Optimization**: Nelder-Mead simplex and Genetic Algorithm optimizers
- **Material Database**: AS4/3501-6, T300/5208, IM7/8552, E-glass/Epoxy, Kevlar/Epoxy

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

## Getting Started

```bash
flutter pub get
flutter run -d linux   # or -d chrome for web
```

## Running Tests

```bash
flutter test
flutter test test/core/math/
flutter test test/integration/
```

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
(-ω²[M] + [K](1+ig) - q·[Q]){η} = 0
```
Sweep velocity V; flutter = minimum V where Im(p)=0 and Re(p) gives flutter frequency.

### Divergence Condition

Static aeroelastic instability when:
```
det([K] - q_D·[A_s]) = 0
q_D = λ_min([K]^{-1}·[A_s])
V_D = √(2·q_D/ρ)
```

## License

MIT License - See LICENSE file for details.
