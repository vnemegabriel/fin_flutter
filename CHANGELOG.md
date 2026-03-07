# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] — 2026-03-07

### Added

**Core Solvers**
- Pure-Dart dense matrix library (`Matrix`, `ComplexMatrix`, `SparseMatrix`) with LU decomposition, inverse, eigenvalue iteration
- Classical Lamination Theory: full [A][B][D] matrix computation for arbitrary angle-ply stacks
- Orthotropic ply model with reduced stiffness [Q] and angle-transformation to [Q̄]
- Failure criteria: Tsai-Wu and Tsai-Hill with reserve factor computation
- FEA engine with Kirchhoff DKQ and Mindlin MITC4 quadrilateral plate elements
- Structured trapezoidal fin mesh generator
- Clamped-root boundary conditions via penalty method (K_ii += 10³⁰)
- Generalized eigenvalue solver (inverse iteration with M-orthogonal deflation + Rayleigh quotient)
- Vortex Lattice Method (VLM): horseshoe vortex, Biot-Savart law, AIC assembly, Kutta-Joukowski lift
- Prandtl-Glauert compressibility correction for subsonic flow (M < 0.85)
- Aeroelastic coupling: displacement interpolation matrix H, modal projection Q_modal = ΦᵀHᵀAICHΦ
- Flutter solver (U-g method): velocity sweep, damping indicator g_i, zero-crossing detection
- Divergence solver: minimum positive eigenvalue of [K]⁻¹[A_s] via inverse power iteration

**Optimization**
- Nelder-Mead simplex optimizer (sync* generator for UI streaming)
- Genetic Algorithm with SBX crossover and polynomial mutation
- Multi-objective loss function with penalty constraints

**Integration**
- OpenRocket (.ork) file parser: ZIP extraction, XML parsing, trapezoidfinset geometry + flight data
- ISA 1976 Standard Atmosphere model (sea level to 30 km)

**Material Database**
- AS4/3501-6 carbon/epoxy
- T300/5208 carbon/epoxy
- IM7/8552 carbon/epoxy
- E-glass/Epoxy
- Kevlar/Epoxy

**Flutter UI**
- Geometry screen with live planform preview
- Materials screen: ply selection, angle slider, CLT display
- Flight condition screen with ISA auto-computation
- Analysis configuration screen (mesh density, element type)
- Results screen: V-g diagram, mode shapes, flutter/divergence speed cards, stability margin status
- Optimization screen: real-time convergence chart
- OpenRocket import workflow

**Documentation**
- Theory reference (`docs/theory.md`): 12-section mathematical derivation
- Benchmark test cases (`docs/test_cases.md`): 9 cases with analytical solutions
- User guide (`docs/user_guide.md`): screen-by-screen walkthrough
- Architecture guide (`docs/architecture.md`): module dependency graph, data-flow table, design decisions
- API reference (`docs/api_reference.md`): headless use without UI
- Contributing guide (`docs/contributing.md`): extension points, code style
- Installation guide (`docs/INSTALL.md`): Linux, macOS, Windows, Web

[1.0.0]: https://github.com/vnemegabriel/fin_flutter/releases/tag/v1.0.0
