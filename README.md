# fin_flutter

An aerospace engineering tool for the analysis of aerodynamic **flutter** and **divergence** of rocket fins — structural aeroelastic instability phenomena.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++](https://img.shields.io/badge/C%2B%2B-17-blue)](https://en.cppreference.com/w/cpp/17)
[![Python](https://img.shields.io/badge/Python-3.10+-blue)](https://www.python.org)
[![CMake](https://img.shields.io/badge/CMake-3.16+-blue)](https://cmake.org)

---

## Overview

**fin_flutter** is a computational platform for aeroelastic analysis of rocket fin systems. It combines:

- **Structural Dynamics** (FEA): Kirchhoff DKQ and Mindlin MITC4 plate elements
- **Aerodynamics** (CFD): Vortex Lattice Method (VLM) with Prandtl-Glauert compressibility
- **Composite Materials**: Classical Lamination Theory (CLT), full [A][B][D] matrix
- **Aeroelastic Stability**: U-g flutter method, static divergence eigenvalue analysis
- **Optimization**: Nelder-Mead simplex, genetic algorithm design space exploration
- **Import/Export**: OpenRocket `.ork` file parser for rocket fin geometry

The **C++ core** provides high-performance numerical solvers. The **Python layer** orchestrates analysis, plotting, and CLI tooling.

---

## Quick Start

### Prerequisites

- **C++17** compiler (g++, clang, MSVC)
- **CMake** ≥ 3.16
- **Eigen3** ≥ 3.3
- **Python** ≥ 3.10 (for orchestration layer)

### Build C++ Core

```bash
cd core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run benchmark tests
./build/tests/run_tests
```

### Setup Python Tools (coming soon)

```bash
cd core/python
pip install -e ".[dev]"
pytest                           # run test suite
python -m fin_flutter.cli --help # CLI entry point
```

See [INSTALL.md](docs/INSTALL.md) for detailed platform-specific instructions.

---

## Features

### Computational Domains

| Domain | Method | Status |
|--------|--------|--------|
| **FEA** | Kirchhoff DKQ, Mindlin MITC4 plate elements | ✓ Implemented |
| **VLM** | Horseshoe vortices, AIC matrix, Kutta-Joukowski | ✓ Core working, sign convention fixed |
| **CLT** | [A][B][D] matrices, ply transformations | ✓ Passing benchmark tests |
| **ISA** | 1976 Standard Atmosphere model | ✓ Passing benchmark tests |
| **Aeroelastic Coupling** | Modal projection, U-g flutter, divergence | ✓ Framework in place |
| **Optimization** | Nelder-Mead, GA design search | ⚠ Framework ready, not yet integrated |

### Material Database

- **AS4/3501-6**: Carbon/epoxy (aerospace standard)
- **T300/5208**: Carbon/epoxy (legacy)
- **IM7/8552**: Carbon/epoxy (high-modulus)
- **E-glass/Epoxy**: Fiberglass (cost-effective)

---

## Repository Structure

```
fin_flutter/
├── core/cpp/                    — C++ computational library
│   ├── include/fin_flutter/
│   │   ├── math/                — Matrix algebra, Eigen types
│   │   ├── materials/           — CLT, orthotropic plies
│   │   ├── fea/                 — Eigenvalue solvers
│   │   ├── cfd/                 — VLM, Biot-Savart, AIC
│   │   ├── aeroelastic/         — Coupling, flutter, divergence
│   │   └── models/              — Geometry, flight conditions
│   ├── tests/                   — Benchmark suite (23 cases)
│   └── CMakeLists.txt
│
├── core/python/                 — Python orchestration (in progress)
│   ├── fin_flutter/
│   │   ├── pipeline.py          — Analysis orchestrator
│   │   ├── plotting.py          — V-g diagrams, mode visualization
│   │   ├── openrocket.py        — .ork file parser
│   │   └── cli.py               — Command-line interface
│   ├── tests/                   — pytest suite
│   └── pyproject.toml
│
├── docs/                        — Technical documentation
│   ├── theory.md                — Mathematical reference
│   ├── test_cases.md            — Benchmark cases with solutions
│   ├── api_reference.md         — C++/Python API docs
│   ├── architecture.md          — Design and data-flow
│   └── INSTALL.md               — Build instructions
│
├── assets/                      — Material databases, tables
├── Recursos/                    — Reference papers (PDFs)
│
├── CLAUDE.md                    — AI assistant session instructions
├── TASK.md                      — Current development status
├── LICENSE                      — MIT License
└── README.md                    — This file
```

---

## Test Status

**Benchmark Suite** (23 test cases):

| Domain | Tests | Status |
|--------|-------|--------|
| Classical Lamination Theory (CLT) | 6 | ✓ All pass |
| ISA 1976 Atmosphere | 9 | ✓ All pass |
| Vortex Lattice Method (VLM) | 1 | ⚠ Sign fixed, magnitude TBD |

**Known Issues:**

- **VLM Case 4** (flat plate, AR=5, α=5°): CL sign is now correct but magnitude undersized by ~50×. Under investigation.

Run tests:

```bash
cd core/cpp
./build/tests/run_tests
```

---

## API Example

### C++ (Header-Only)

```cpp
#include <fin_flutter/materials/clt_calculator.hpp>
#include <fin_flutter/models/flight_condition.hpp>
#include <fin_flutter/cfd/vortex_lattice.hpp>

using namespace ff;

// 1. Define laminate
std::vector<LaminatePly> plies = {
  {OrthotropicPly::as4_3501_6(), 0},
  {OrthotropicPly::as4_3501_6(), 45},
  {OrthotropicPly::as4_3501_6(), -45},
  {OrthotropicPly::as4_3501_6(), 90}
};

// 2. Compute CLT matrices
CLTCalculator clt;
auto laminate = clt.compute(plies);
std::cout << "A11 = " << laminate.A(0,0) / 1e6 << " MN/m\n";

// 3. Define flight condition
auto fc = FlightCondition::from_isa(altitude_m, velocity_ms);

// 4. Run VLM
FinGeometry fin{0.5, 0.2, 0.2, 0.0, 0.0, 1};
VortexLattice vlm;
auto result = vlm.solve(fin, fc, 8, 12, alpha_rad);
std::cout << "CL = " << result.CL << "\n";
```

### Python (Coming Soon)

```python
from fin_flutter import CLTCalculator, FlightCondition, VortexLattice

# CLT analysis
laminate = CLTCalculator.from_stack([
    ('AS4/3501-6', 0),
    ('AS4/3501-6', 45),
    ('AS4/3501-6', -45),
    ('AS4/3501-6', 90),
])
print(f"A11 = {laminate.A[0,0]/1e6:.1f} MN/m")

# VLM analysis
fc = FlightCondition.from_isa(altitude_m=0, velocity_ms=50)
fin = FinGeometry(span=0.5, root_chord=0.2, tip_chord=0.2)
vlm = VortexLattice(fin, fc)
result = vlm.solve(alpha_deg=5, chordwise_panels=8, spanwise_panels=12)
print(f"CL = {result.CL:.3f}")
```

---

## Documentation

| Document | Content |
|----------|---------|
| [docs/theory.md](docs/theory.md) | Complete mathematical derivations: CLT, FEA, VLM, aeroelastic coupling, flutter/divergence methods |
| [docs/test_cases.md](docs/test_cases.md) | 9 benchmark cases with analytical reference solutions and test tolerances |
| [docs/architecture.md](docs/architecture.md) | Design patterns, module dependency graph, data-flow through the pipeline |
| [docs/api_reference.md](docs/api_reference.md) | C++ Doxygen API reference and Python function signatures |
| [docs/INSTALL.md](docs/INSTALL.md) | Build from source on Linux, macOS, Windows; dependency management |
| [CLAUDE.md](CLAUDE.md) | Session instructions for AI-assisted development (token efficiency, coding standards) |

---

## Development Status

**Current Focus**: Fix VLM Case 4 magnitude undersizing (CL = 0.0077 vs expected 0.36–0.42)

**Completed**:
- ✓ C++ core library structure
- ✓ Biot-Savart vortex-induced velocity (corrected formula)
- ✓ Horseshoe vortex sign convention (trailing vortex Γ fix)
- ✓ CLT laminate analysis with full [A][B][D] matrices
- ✓ ISA 1976 atmosphere model (6 layers)
- ✓ FEA eigenvalue solver framework
- ✓ Aeroelastic coupler and modal projection
- ✓ Flutter U-g and divergence solver frameworks
- ✓ CMake build system with Eigen integration
- ✓ Benchmark test suite (23 cases)

**In Progress**:
- Python orchestration layer (pipeline.py, cli.py, plotting.py)
- VLM magnitude scaling investigation
- Integration test (end-to-end pipeline)

**Future**:
- OpenRocket `.ork` file parser
- Design optimization (Nelder-Mead, genetic algorithm)
- Visualization and plotting (V-g diagrams, mode shapes)
- Web interface

---

## Contributing

See [docs/contributing.md](docs/contributing.md) for guidelines on:
- Adding new materials to the database
- Implementing FEA elements
- Extending the VLM with corrections (sweep, dihedral)
- Adding benchmark test cases

For development notes, see [CLAUDE.md](CLAUDE.md) for session structure and coding standards.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

Copyright © 2025 fin_flutter contributors

---

## References

**Key Papers & Textbooks:**

- Katz, J. & Plotkin, A. (2001). *Low-Speed Aerodynamics* (2nd ed.). Cambridge University Press.
- Bisplinghoff, R.L., Ashley, H. & Halfman, R.L. (1955). *Aeroelasticity*. Addison-Wesley.
- Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells* (2nd ed.). CRC Press.
- NACA Report 1235: The 1976 Standard Atmosphere (U.S. Standard Atmosphere).

See [Recursos/](Recursos/) for full reference papers and datasheets.
