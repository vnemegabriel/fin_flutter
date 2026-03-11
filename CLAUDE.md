# CLAUDE.md — Session Instructions for fin_flutter

This file is read automatically at the start of every Claude Code session.
It establishes non-negotiable architectural rules and efficient working patterns.

---

## Token Efficiency — Mandatory

- **Be concise in all responses.** No filler phrases, no restating the task, no summaries of what was just done.
- **Lead with action or answer.** Reasoning only when it changes the decision.
- **Do not read files you don't need.** Load only what the task requires (see Context Loading Guide).
- **Do not output code you haven't changed.** Show diffs or targeted snippets, not full file rewrites.
- **One tool call per logical unit.** Batch parallel reads; do not chain unnecessary sequential reads.
- **Skip confirmations for reversible local actions.** Just do it and report the result.

---

## Repository Overview

**fin_flutter** is an aerospace engineering tool for the analysis of aerodynamic **flutter** and **divergence** of rocket fins — structural aeroelastic instability phenomena.

**Computational domains:**

| Domain | Method |
|--------|--------|
| Structural dynamics | FEA — Kirchhoff DKQ and Mindlin MITC4 plate elements |
| Aerodynamics | Vortex Lattice Method (VLM) with Prandtl-Glauert compressibility |
| Composite materials | Classical Lamination Theory (CLT), full [A][B][D] matrix |
| Aeroelastic stability | U-g flutter method, static divergence eigenvalue |
| Design optimization | Nelder-Mead simplex, genetic algorithm |
| Geometry import | OpenRocket `.ork` file parser |

**Language strategy (in priority order):**

| Priority | Language | Role |
|----------|----------|------|
| 1 | **C++** | Performance-critical numerical core: solvers, matrix ops, FEA, VLM, CLT |
| 2 | **Python** | Orchestration, scripting, plotting, CLI tooling, testing |
| 3 | Julia | Numerical experimentation, rapid prototyping |
| 4 | MATLAB | Legacy reference implementations only |

---

## Hard Rules — Never Violate Without Explicit Instruction

### Numerical Correctness

- **Do not change solver algorithms** without referencing `docs/theory.md` and validating against `docs/test_cases.md`.
- **All numerical changes must pass benchmark validation.** Compare output against analytical solutions in `docs/test_cases.md` before committing.
- **Physical constants, material properties, and CLT formulation are constraints.** Do not alter without explicit instruction.

### Language Rules

- **New code goes in C++ or Python.** Do not write solvers or analysis logic in any other language unless explicitly requested.
- **C++ for computation; Python for orchestration and tooling.** Do not invert this.
- **Julia and MATLAB are offline/experimental only.** Never in production or CI paths.

---

## Engineering Programming Standards

- **Modular, single-responsibility code.** One class/module per file. Split when a file exceeds ~200 lines.
- **Documented numerical methods.** Every non-trivial formula must cite the source equation and reference:
  - C++: `// Eq. 4.12 — Theodorsen U-g method [Bisplinghoff, Ashley & Halfman, 1955]`
  - Python: `# Eq. 4.12 — Theodorsen U-g method [Bisplinghoff, Ashley & Halfman, 1955]`
- **Concise implementations.** Clarity over abstraction. Three explicit lines beat a premature helper.
- **Test against analytical solutions.** Every new numerical routine requires a test with a known closed-form result from `docs/test_cases.md`.

### C++ Standards

- C++17 minimum. Use **Eigen** for dense/sparse linear algebra.
- Headers in `include/`, implementations in `src/`. One class per header.
- **CMake** build system. No hand-written Makefiles.
- **Doxygen** comments on all public interfaces (`@brief`, `@param`, `@return`, equation reference).

### Python Standards

- Python 3.10+. Use **NumPy/SciPy** for numerics, **Matplotlib** for plotting.
- **Type annotations** required on all public functions.
- **Google-style docstrings** with parameter types, return values, and equation references.
- `pyproject.toml` for packaging; **pytest** for tests.

---

## Documentation Standards

- Every public C++ interface (`include/`) requires a Doxygen comment with `@brief`, `@param`, `@return`, and equation reference.
- Every public Python function requires a Google-style docstring with types and equation reference.
- When adding a new analysis feature, update `docs/api_reference.md` and `docs/theory.md`.
- New benchmark cases must be added to `docs/test_cases.md` with the analytical reference and full citation.

---

## Architectural Map

```
core/cpp/                    — C++ computational library
  include/                   — Public headers (Doxygen-documented)
    math/                    — Matrix, complex arithmetic, sparse solvers
    materials/               — CLT, orthotropic ply, failure criteria
    fea/                     — Mesh, elements, assembly, eigensolvers
    cfd/                     — VLM, AIC matrix, horseshoe vortex, compressibility
    aeroelastic/             — Aeroelastic coupler, flutter U-g, divergence
    optimization/            — Nelder-Mead, genetic algorithm
  src/                       — Implementations
  tests/                     — Unit + benchmark tests (Google Test or Catch2)
  CMakeLists.txt
core/python/                 — Python orchestration and tooling
  fin_flutter/               — Package source
    pipeline.py              — Analysis pipeline orchestrator
    plotting.py              — V-g diagrams, mode shapes
    openrocket.py            — .ork file parser
    cli.py                   — Command-line interface
  tests/                     — pytest suite
  pyproject.toml
docs/                        — theory.md, architecture.md, api_reference.md, test_cases.md, INSTALL.md
assets/                      — Material database, ISA atmosphere tables, example projects
Recursos/                    — Reference papers (PDFs)
```

---

## Context Loading Guide

| Task type | Load these files |
|-----------|-----------------|
| FEA / structural solver | `core/cpp/include/fea/`, `core/cpp/src/fea/`, `docs/theory.md` (FEA sections) |
| Aerodynamics / VLM | `core/cpp/include/cfd/`, `core/cpp/src/cfd/`, `docs/theory.md` (VLM sections) |
| Aeroelastic coupling / flutter | `core/cpp/include/aeroelastic/`, `docs/theory.md` (flutter sections) |
| Materials / CLT | `core/cpp/include/materials/`, `docs/theory.md` (CLT sections) |
| Matrix / linear algebra | `core/cpp/include/math/` |
| Optimization | `core/cpp/include/optimization/` or `core/python/fin_flutter/` |
| OpenRocket import | `core/python/fin_flutter/openrocket.py` |
| CLI / scripting | `core/python/fin_flutter/cli.py`, `core/python/fin_flutter/pipeline.py` |
| Benchmarks / validation | `docs/test_cases.md`, `core/cpp/tests/`, `core/python/tests/` |
| Architecture overview | `docs/architecture.md` |

Do not load all files for every session. Load the minimum set.

---

## Task Scoping Rules

- **One domain per session where possible.** Structural, aerodynamic, and aeroelastic work should be separate sessions unless the task explicitly spans multiple domains.
- **Write a TASK.md before starting** anything that touches 3 or more files. Include: current state, required changes, acceptance criteria, out-of-scope items.
- **State what NOT to change.** Explicitly declare out-of-scope areas to prevent scope drift.
- **Batch related small changes.** Do not batch unrelated work.
- **Validate numerics before any commit.** Changes to numerical code require benchmark validation against `docs/test_cases.md`.

---

## Building and Running

```bash
# C++ core
cd core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/tests/run_tests          # run benchmark suite

# Python tools
cd core/python
pip install -e ".[dev]"
pytest                           # full test suite
python -m fin_flutter.cli --help # CLI entry point
```
