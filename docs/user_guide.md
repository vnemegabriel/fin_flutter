# User Guide

This guide walks through every screen in Fin Flutter Analyzer, explaining inputs,
outputs, and tips for getting accurate results.

---

## 1. Requirements and Launch

**Install:** See [INSTALL.md](INSTALL.md) for platform setup.

```bash
flutter pub get
flutter run -d linux      # or -d chrome, -d macos, -d windows
```

The app opens on the **Home** screen.

---

## 2. Home Screen

Displays a list of saved projects. Each project stores geometry, material, flight
condition, and the most recent analysis result.

**Actions:**
- **New Project** — opens the Geometry screen with blank defaults.
- **Open Project** — loads a previously saved project.
- **Import .ork** — launches a file picker to load an OpenRocket file (see §10).

---

## 3. Geometry Screen

Define the trapezoidal fin planform.

| Parameter | Units | Typical Range |
|-----------|-------|---------------|
| Span | m | 0.05 – 0.5 |
| Root chord | m | 0.05 – 0.3 |
| Tip chord | m | 0 – root chord |
| Sweep length | m | 0 – span |
| Thickness | m | 0.001 – 0.01 |

The **planform preview** updates in real time as you drag sliders.

**Derived properties** (shown below the preview):
- Planform area, aspect ratio, sweep angle, taper ratio, mean aerodynamic chord.

**Tips:**
- Keep thickness ≤ 5% of chord length for thin-plate theory to be valid.
- A taper ratio of 0.4–0.6 is typical for amateur rockets.
- Very low aspect ratio (AR < 1) fins produce less accurate VLM results.

---

## 4. Materials Screen

### 4.1 Ply Material Selection

Choose a prepreg material from the database:

| Material | E₁ (GPa) | E₂ (GPa) | G₁₂ (GPa) |
|----------|----------|----------|-----------|
| AS4/3501-6 | 142 | 10.3 | 7.2 |
| T300/5208 | 132 | 10.8 | 5.65 |
| IM7/8552 | 171 | 9.08 | 5.29 |
| E-glass/Epoxy | 45 | 12.0 | 5.5 |
| Kevlar/Epoxy | 76 | 5.5 | 2.3 |

### 4.2 Laminate Stack

Add plies using the **+ Add Ply** button. For each ply, set:
- **Angle (°):** fiber orientation relative to spanwise direction. 0° = spanwise.
- **Thickness:** usually set automatically from the material's nominal ply thickness.

**Preset layups:**
- **Unidirectional [0]ₙ** — maximum bending stiffness (minimum flutter speed for given thickness).
- **Cross-ply [0/90]ₙs** — intermediate bending and twisting stiffness.
- **Quasi-isotropic [0/45/−45/90]ₛ** — balanced in-plane response; recommended starting point.

### 4.3 CLT Display

After building the stack, the computed [A], [B], [D] matrices and derived engineering
constants (E_x, E_y, ν_xy, G_xy) are displayed.

**Design guidance:**
- Aim for B = 0 (symmetric laminate) to avoid bending-extension coupling.
- D₁₁ (spanwise bending stiffness) is the most important parameter for flutter resistance.
- Increase the number of 0° plies to raise D₁₁; add ±45° plies to raise D₆₆ and reduce flutter susceptibility.

---

## 5. Flight Conditions Screen

| Parameter | Description |
|-----------|-------------|
| Altitude | Flight altitude in metres (ISA model) |
| Mach number | Flight Mach number |

ISA values (density ρ, speed of sound a, dynamic pressure q) are computed automatically
and displayed below the inputs.

**Warning:** A banner appears if M ≥ 0.85. The Prandtl-Glauert correction used by the VLM
becomes inaccurate above M = 0.85 (transonic regime). For transonic analysis, use a
dedicated Euler/Navier-Stokes solver.

**Max flight velocity** — used to compute stability margins. Set to the maximum velocity
the fin will experience in flight (including motor burnout and coast phase peak velocity).

---

## 6. Analysis Configuration Screen

### 6.1 Mesh Density

| Mesh | Spanwise | Chordwise | Total DOF (approx.) | Eigensolve time |
|------|----------|-----------|---------------------|-----------------|
| Coarse (fast) | 4 | 3 | ~75 | < 0.2 s |
| Medium | 6 | 4 | ~175 | ~0.5 s |
| Fine | 12 | 8 | ~700 | ~3–5 s |
| Very fine | 20 | 12 | ~1800 | ~15–30 s |

For initial design, use **Medium**. For final validation, use **Fine**.

### 6.2 Element Type

- **Kirchhoff DKQ** (default) — thin plate; valid when h/c < 0.05. Faster.
- **Mindlin MITC4** — thick plate; use when h/c ≥ 0.05 or for higher-order mode accuracy.

### 6.3 Number of Modes

More modes capture higher-frequency flutter coupling. 2–4 modes is usually sufficient;
adding more modes increases analysis time without changing V_flutter significantly.

### 6.4 VLM Panels

Default: 4 chordwise × 6 spanwise = 24 panels. Increasing to 8×12 improves CL accuracy
by ~2% but barely changes flutter/divergence speeds (the aeroelastic coupling is dominated
by the lowest aerodynamic mode).

---

## 7. Running the Analysis

Press **Run Analysis**. The compute service dispatches to a background isolate.

Progress indicator shows: FEA → VLM → Coupling → Flutter → Divergence.

Results appear on the **Results** screen when complete (typically 0.5–5 s).

---

## 8. Results Screen

### 8.1 Summary Cards

| Card | Value |
|------|-------|
| Flutter Speed | V_flutter (m/s), flutter margin (V_F/V_max) |
| Divergence Speed | V_divergence (m/s), divergence margin (V_D/V_max) |
| Flutter Frequency | Frequency of the fluttering mode at V_F (Hz) |
| Safety Status | SAFE / WARNING / MARGINAL / CRITICAL (colour coded) |

### 8.2 V-g Diagram

The V-g diagram plots the damping indicator g_i(V) for each structural mode against velocity.

- **Horizontal axis:** velocity (m/s)
- **Vertical axis:** damping g (dimensionless)
- **g = 0 line:** stability boundary. The velocity where g first crosses zero is V_flutter.
- Each mode is a separate coloured curve.

**Interpretation:**
- Modes whose g never reaches 0 are stable throughout the velocity range.
- The mode that crosses g = 0 first determines V_flutter.
- If two modes coalesce (similar frequencies), coupled flutter may occur — visible as two curves merging.

### 8.3 Natural Frequencies Tab

Lists the computed structural natural frequencies (Hz) and mode shapes. Mode 1 is typically
the first bending mode; Mode 2 is the first torsional mode.

### 8.4 Mode Shape Tab

An interactive 3D view of each mode shape deformation. The magnitude scale is normalised.

---

## 9. Optimization Screen

### 9.1 Algorithm Selection

| Algorithm | Best for | Notes |
|-----------|----------|-------|
| Nelder-Mead | Few variables (< 10), smooth objective | Fast; deterministic |
| Genetic Algorithm | Many variables, multimodal objective | Stochastic; slower but more global |

### 9.2 Design Variables

Check the parameters to optimize:
- Span, root/tip chord, sweep — geometric
- Ply angles and thicknesses — material

Narrower bounds lead to faster convergence but may miss a global optimum.

### 9.3 Convergence Chart

The chart streams each iteration's best objective value in real time.

**Nelder-Mead:** convergence is typically monotonic; stop criteria is standard deviation of
vertex values < 10⁻⁶.

**GA:** expect noisy early generations; the curve smooths as the population converges.
If the objective plateaus far from zero, widen the design variable bounds or increase
population size.

### 9.4 Applying the Optimum

Press **Apply Best** to populate the Geometry and Materials screens with the optimized
values, then re-run the analysis to confirm.

---

## 10. OpenRocket Import

1. In any screen, press **File → Open .ork** (or tap the import button on the Home screen).
2. Select a `.ork` file from your file system.
3. The parser extracts the first `trapezoidfinset` found in the XML.
4. Geometry screen is auto-populated with span, root chord, tip chord, sweep, and thickness.
5. Flight data (if present) sets the altitude and estimated max velocity.
6. Review the imported values and adjust if needed (e.g., material is not imported).

**Supported OpenRocket versions:** 1.5 – 23.09.

---

## 11. Exporting Results

Press **Export → PDF Report** from the Results screen. The PDF includes:
- Fin geometry diagram
- Material and laminate summary
- Flight condition
- Flutter and divergence speeds with stability margin status
- V-g diagram image

---

## 12. Interpretation Guidance

### Is my design safe?

- Target a flutter margin ≥ 1.5 (V_flutter ≥ 1.5 × V_max) for amateur rockets.
- The divergence margin should also be ≥ 1.5; in practice V_D > V_F for most fin designs.
- If either margin is CRITICAL (< 1.0), the fin will flutter or diverge before reaching max velocity.

### Conservative vs. optimistic settings

| Setting | Conservative | Optimistic |
|---------|-------------|------------|
| Element type | Mindlin MITC4 | Kirchhoff DKQ |
| Mesh | Fine | Coarse |
| VLM panels | 8×12 | 4×4 |
| Modes | 4 | 2 |

For a final design review, always use conservative settings.

### Common issues

| Symptom | Likely cause |
|---------|-------------|
| V_flutter very low (< 50 m/s) for CFRP fin | Laminate is all-0° (no torsional stiffness) |
| V_divergence returns "No divergence found" | A_s has no positive eigenvalue (divergence-stable) |
| Flutter frequency = 0 | Mode with zero frequency found; check BC application |
| NaN in results | Degenerate mesh (zero-length edge) or zero-thickness laminate |
