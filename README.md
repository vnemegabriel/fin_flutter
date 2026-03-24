# FalconLAUNCH VI Fin Flutter Analysis

A comprehensive aeroelastic analysis toolkit for predicting and validating fin flutter behavior in supersonic rocket applications. This project combines classical laminate theory with empirical flutter prediction methods to assess structural stability across the vehicle envelope.

## Overview

Fin flutter is a critical aeroelastic phenomenon in high-speed flight where aerodynamic loads induce self-excited structural oscillations. This analysis framework predicts the flutter speed boundary and verifies structural margins against max-q conditions using validated methods from NACA and aerospace industry standards.

### Key Capabilities

- **Classical Laminate Theory (CLPT)**: Analyze composite fin layups with arbitrary fiber orientations
- **Aeroelastic Tailoring**: Apply global beta-rotation (Weisshaar 1981) for optimized stiffness tuning
- **Subsonic & Supersonic Flutter Prediction**: Empirical correlations spanning Mach 0.3–3.0+
- **Sweep Corrections**: Account for leading-edge sweep using established correction factors
- **Flight Envelope Integration**: Assess flutter safety factors across altitude and Mach profiles
- **Composite Micromechanics**: Halpin-Tsai modeling for effective ply properties

## Theory & References

This toolkit is built on established aeroelastic theory and composite mechanics:

### Primary References

1. **Martin, H.C.** (1958) *NACA TN 4197* — Classical fin flutter boundary estimation
2. **Ackeret, J.** (1925) *NACA TM 317* — Prandtl-Ackeret supersonic lift coefficients
3. **Bisplinghoff, R.L., Ashley, H. & Halfman, R.L.** (1955) *Aeroelasticity* — Sweep correction methodology (Section 5.5)
4. **Weisshaar, T.A.** (1981) *J. Aircraft* 18(8):669–676 — Aeroelastic tailoring via global rotation
5. **Halpin, J.C. & Tsai, S.W.** (1969) *AFML-TR-67-423* — Composite micromechanics
6. **Jones, R.M.** (1999) *Mechanics of Composite Materials* (2nd ed.) — Classical laminate theory
7. **Shirk, M.H., Hertz, T.J. & Weisshaar, T.A.** (1986) *J. Aircraft* 23(1):6–21 — Composite structural tailoring
8. **AIAA S-080** (1999) & **MIL-A-8870C** (1993) — Industry standards

## Installation

### Requirements

- Python 3.8+
- Standard library dependencies: `math`, `json`, `argparse`, `pathlib`

### Setup

```bash
# Clone the repository
git clone <repository-url>
cd fin_flutter

# No external package dependencies required
python3 flutterEstimate_v3.py
```

## Quick Start

### 1. Laminate Analysis: In-Plane Stiffness

Analyze composite fin layups and generate laminate properties:

```bash
# Default analysis: AR1 architecture, fiber volume fraction 50%, full beta sweep
python3 inplaneG_v5.py

# Single beta value (e.g., 20° orientation)
python3 inplaneG_v5.py --beta 20

# Custom fiber volume fraction
python3 inplaneG_v5.py --vf 0.55 --beta 20

# Scale to target thickness (6 mm total)
python3 inplaneG_v5.py --thickness 6.0

# Export to JSON for flutter pipeline
python3 inplaneG_v5.py --beta 20 --json laminate.json

# Silent JSON export (suppress console output)
python3 inplaneG_v5.py --quiet --json laminate.json
```

### 2. Flutter Boundary Prediction

Predict flutter speed and verify safety margins:

```bash
# Default analysis (Mach 1.942 at 1462 m altitude)
python3 flutterEstimate_v3.py

# Use laminate properties from JSON
python3 flutterEstimate_v3.py --json laminate.json

# Manual in-plane shear stiffness (D66 [N⋅m])
python3 flutterEstimate_v3.py --d66 205.96

# Custom flight conditions
python3 flutterEstimate_v3.py --altitude 1462 --mach 1.942

# Sweep table across altitudes
python3 flutterEstimate_v3.py --sweep-table

# Flight envelope analysis from CSV
python3 flutterEstimate_v3.py --csv flight_data.csv
```

### 3. Integrated Workflow

Complete pipeline from laminate design to flutter validation:

```bash
# Step 1: Design laminate at optimal fiber angle
python3 inplaneG_v5.py --beta 25 --thickness 7.5 --json fin_layup.json

# Step 2: Predict flutter speed with designed properties
python3 flutterEstimate_v3.py --json fin_layup.json --sweep-table

# Step 3: Assess against flight envelope
python3 flutterEstimate_v3.py --json fin_layup.json --csv flight_data.csv
```

## Key Parameters

### `inplaneG_v5.py` (Laminate Analysis)

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `--vf` | 0.50 | — | Fiber volume fraction (0.0–0.65) |
| `--beta` | *sweep* | deg | Off-axis fiber angle |
| `--layup` | ar1 | — | Half-stack architecture (ar1, more_db, more_ga, equal) |
| `--thickness` | *nominal* | mm | Scale laminate to target total thickness |
| `--json` | — | path | Write laminate properties to JSON |
| `--quiet` | false | — | Suppress console output |

**Layup Architectures:**
- **ar1**: Aspect-ratio-driven (baseline)
- **more_db**: Higher D22 (in-plane transverse stiffness)
- **more_ga**: Higher G12 (in-plane shear stiffness)
- **equal**: Symmetric in-plane properties

### `flutterEstimate_v3.py` (Flutter Analysis)

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `--d66` | *calculated* | N⋅m | In-plane shear stiffness override |
| `--altitude` | 1462 | m | Analysis altitude (max-q condition) |
| `--velocity` | *ISA-derived* | m/s | Rocket velocity |
| `--mach` | 1.942 | — | Mach number |
| `--span` | 0.160 | m | Fin semi-span |
| `--cr` | 0.300 | m | Root chord |
| `--ct` | 0.150 | m | Tip chord |
| `--sweep` | 57.4 | deg | Leading-edge sweep angle |
| `--sweep-off` | false | — | Disable sweep correction |
| `--super-off` | false | — | Disable supersonic correction |
| `--fsf-req` | 1.50 | — | Required flutter safety factor |
| `--csv` | — | path | Flight data CSV (overrides altitude/velocity) |

## Project Structure

```
fin_flutter/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── flutterEstimate_v3.py             # Flutter speed prediction (subsonic/supersonic)
├── inplaneG_v5.py                    # Laminate analysis & CLT
├── flight_data.csv                   # Example flight envelope data
├── layup_viewer_v5.html              # Interactive laminate visualization
├── core/
│   └── python/
│       ├── fin_flutter/
│       │   └── flight_data.py        # Flight envelope utilities
│       └── tests/
│           └── test_flutter_equations.py  # Unit tests
└── Recursos/                          # Reference materials (PDFs)
    ├── Aeroelastic Tailoring — Theory, Practice, and Promise.pdf
    ├── Libro-Dowell-A Modern Course in Aeroelasticity.pdf
    ├── Libro-Bisplinghoff-Aeroelasticity-1983.pdf
    ├── Introduction to Composite Materials Design.pdf
    ├── Mechanicsofcompositematerials-Jones1999.pdf
    └── [Additional technical references]
```

## Core Modules

### `fin_flutter.flight_data`

Utilities for reading and processing flight telemetry:

```python
from fin_flutter.flight_data import read_flight_data

# Load flight data CSV
data = read_flight_data("flight_data.csv")

# Identify max-velocity condition
max_v_row = data[data['velocity'].idxmax()]
altitude_max_q = max_v_row['altitude']
velocity_max_q = max_v_row['velocity']
```

## Running Tests

```bash
python3 -m pytest core/python/tests/ -v
```

## Interpretation & Outputs

### Flutter Speed (V_f)

The velocity at which aeroelastic instability occurs. Output includes:

- **Subsonic prediction** (M < 0.9): Martin 1958 correlation
- **Transonic correction** (0.9 < M < 1.3): Interpolated
- **Supersonic prediction** (M > 1.3): Ackeret-corrected formula
- **Sweep-corrected flutter speed**: Accounts for leading-edge sweep reduction in effective stiffness

### Safety Factor (FSF)

$$\text{FSF} = \frac{V_f}{\max(V_{\text{flight envelope}})}$$

**Acceptable margins:**
- FSF ≥ 1.50 (typical requirement)
- FSF ≥ 1.15 (minimum acceptable in constrained designs)

### Laminate Output (JSON)

When using `--json`, the output includes:

```json
{
  "layup": "ar1",
  "beta": 25.0,
  "vf": 0.50,
  "thickness_mm": 7.5,
  "properties": {
    "E11": 45230,
    "E22": 12450,
    "G12": 5680,
    "nu12": 0.31,
    "D11": 5428,
    "D22": 1494,
    "D66": 680.4
  }
}
```

## Common Workflows

### Optimize Fiber Angle for Maximum Flutter Speed

```bash
# Sweep beta from 0° to 45°
for b in {0..45..5}; do
  python3 inplaneG_v5.py --beta $b --json fin_$b.json
  python3 flutterEstimate_v3.py --json fin_$b.json
done
# Compare FSF values to identify optimal beta
```

### Validate Against Flight Data

```bash
# Use actual flight profile from telemetry
python3 flutterEstimate_v3.py --csv flight_data.csv --json fin_layup.json
```

### Interactive Layup Design

Open `layup_viewer_v5.html` in a web browser to visualize laminate ply orientations and effective properties interactively.

## Limitations & Assumptions

- **ISA atmosphere**: Valid 0–11 km altitude
- **Linear elasticity**: No material nonlinearity
- **Laminate theory**: Assumes thin shells (thickness << span)
- **Empirical correlations**: Valid for small-aspect-ratio fins (typical rockets)
- **Incompressible approximation**: Transonic regime uses interpolation; not valid near M=1.0 exactly
- **No damping modeling**: Flutter prediction conservative (ignores structural/aerodynamic damping)

## Contributing

For bug reports, feature requests, or technical discussions, please open an issue or submit a pull request.

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

---

**Maintained for FalconLAUNCH VI Project**
Last Updated: 2026-03-24
