# TASK.md — Build C++ / Python Core

## Current State

- Dart/Flutter reference implementation exists in `lib/core/` and `lib/modules/`
- `core/` directory does not exist; no C++ or Python code written
- All domains ported from Dart: math, materials (CLT), CFD (VLM), FEA, aeroelastic, flutter

## Goal

Bootstrap the C++ computational library (`core/cpp/`) and Python tooling (`core/python/`).

## Required Changes

### C++ (`core/cpp/`)
- CMakeLists.txt (root + tests/)
- `include/fin_flutter/math/vec3.hpp` — 3D vector, Biot-Savart operations
- `include/fin_flutter/math/types.hpp` — Eigen matrix typedefs
- `include/fin_flutter/models/fin_geometry.hpp` — trapezoidal planform
- `include/fin_flutter/models/flight_condition.hpp` — ISA 1976, FlightCondition
- `include/fin_flutter/materials/orthotropic_ply.hpp` — Q and Q̄ matrices
- `include/fin_flutter/materials/clt_calculator.hpp` — [A][B][D] matrices
- `include/fin_flutter/cfd/biot_savart.hpp` — Biot-Savart law
- `include/fin_flutter/cfd/vortex_lattice.hpp` — VLM panels, AIC, gamma solve
- `include/fin_flutter/fea/eigenvalue_solver.hpp` — generalized eigenvalue (Eigen)
- `include/fin_flutter/aeroelastic/aeroelastic_coupler.hpp` — modal GAF matrix
- `include/fin_flutter/aeroelastic/flutter_solver.hpp` — U-g sweep
- `include/fin_flutter/pipeline.hpp` — orchestrates full analysis
- `tests/CMakeLists.txt`
- `tests/test_main.cpp` — validates benchmark cases from docs/test_cases.md

### Python (`core/python/`)
- `pyproject.toml`
- `fin_flutter/__init__.py`
- `fin_flutter/pipeline.py` — calls C++ binary, parses JSON output
- `fin_flutter/plotting.py` — V-g diagrams via matplotlib
- `fin_flutter/cli.py` — CLI entry point

## Out of Scope

- Dart/Flutter UI (not modified)
- Optimization module (next session)
- OpenRocket parser (next session)
- pybind11 C++/Python bindings (subprocess JSON interface for now)

## Acceptance Criteria

1. `cmake -B build && cmake --build build` succeeds
2. `./build/tests/run_tests` passes all 3 benchmark cases:
   - CLT [0]₈ AS4/3501-6: A₁₁ = 142.75 MN/m ±0.2%
   - CLT quasi-isotropic: A₁₁ = A₂₂ (isotropy), B = 0
   - ISA: ρ(0) = 1.225 kg/m³, a(0) = 340.29 m/s
3. Python: `python -m fin_flutter.cli --help` runs
