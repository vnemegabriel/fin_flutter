/// @file test_main.cpp
/// @brief Benchmark validation tests for the fin_flutter C++ core.
///
/// Validates against analytical solutions from docs/test_cases.md.
/// Exits 0 on pass, 1 on any failure.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "fin_flutter/materials/clt_calculator.hpp"
#include "fin_flutter/materials/orthotropic_ply.hpp"
#include "fin_flutter/models/flight_condition.hpp"
#include "fin_flutter/cfd/vortex_lattice.hpp"
#include "fin_flutter/models/fin_geometry.hpp"

// ── Test framework (minimal, no dependencies) ─────────────────────────────────

static int g_pass = 0;
static int g_fail = 0;

static void check(bool ok, const char* test_name, const char* detail = "") {
    if (ok) {
        std::printf("  [PASS] %s %s\n", test_name, detail);
        ++g_pass;
    } else {
        std::printf("  [FAIL] %s %s\n", test_name, detail);
        ++g_fail;
    }
}

static void check_rel(double value, double expected, double tol_frac,
                       const char* name) {
    const double err = std::abs(value - expected) / (std::abs(expected) + 1e-30);
    char detail[256];
    std::snprintf(detail, sizeof(detail),
        "(got %.6g, expected %.6g, err=%.4f%%)", value, expected, err*100.0);
    check(err <= tol_frac, name, detail);
}

static void check_abs(double value, double expected, double tol_abs,
                       const char* name) {
    const double err = std::abs(value - expected);
    char detail[256];
    std::snprintf(detail, sizeof(detail),
        "(got %.6g, expected %.6g, abs_err=%.3e)", value, expected, err);
    check(err <= tol_abs, name, detail);
}

// ── Case 1: CLT [0]₈ AS4/3501-6 ─────────────────────────────────────────────
// Reference: docs/test_cases.md §Case 1 — Reddy (2004) Table 3.3.
// Expected: A₁₁ = 142.75 MN/m ±0.2%, B = 0 (±1 mN), D₁₁ = 11.90 N·m ±0.2%.

static void test_clt_unidirectional() {
    std::printf("\nCase 1 — CLT [0]₈ AS4/3501-6\n");

    const ff::OrthotropicPly ply = ff::as4_3501_6(125e-6); // t = 125 µm

    std::vector<ff::LaminatePly> stack;
    for (int i = 0; i < 8; ++i)
        stack.push_back({ply, 0.0}); // all plies at 0°

    ff::CLTCalculator clt;
    const ff::LaminateABD abd = clt.compute(stack);

    // Analytical: Δ = 1 − 0.27 × (0.27×10.3/142) = 0.99471
    // Q₁₁ = 142e9 / 0.99471 = 142.756e9 Pa
    // Q₁₂ = 0.27 × 10.3e9 / 0.99471 = 2.7958e9 Pa
    // A₁₁ = Q₁₁ × h = 142.756e9 × 1e-3 = 142.756e6 N/m
    // A₁₂ = Q₁₂ × h = 2.7958e9 × 1e-3 = 2.7958e6 N/m
    // Note: test_cases.md lists "≈ 2.78 MN/m" which is a rounded approximation.
    check_rel(abd.A(0,0), 142.75e6, 0.002, "A11 [N/m]");
    check_rel(abd.A(0,1),  2.796e6, 0.002, "A12 [N/m]");
    check_rel(abd.A(2,2),   7.20e6, 0.001, "A66 [N/m]");

    // B must be zero for symmetric laminate (rounding ≤ 1e-3 N)
    check_abs(abd.B.norm(), 0.0, 1e-3, "B == 0 [N]");

    // D₁₁ = Q₁₁ × h³/12 = 142.756e9 × (1e-3)³/12 = 11.896 N·m
    check_rel(abd.D(0,0), 11.90, 0.002, "D11 [N·m]");

    // Total thickness
    check_rel(abd.total_thickness, 1e-3, 1e-6, "total_thickness [m]");
}

// ── Case 2: CLT Quasi-isotropic [0/45/−45/90]_s ──────────────────────────────
// Reference: docs/test_cases.md §Case 2 — Reddy (2004) Table 3.7.
// Expected: A₁₁ = A₂₂ (isotropy), A₁₆ ≈ 0, B = 0.

static void test_clt_quasi_isotropic() {
    std::printf("\nCase 2 — CLT [0/45/−45/90]_s quasi-isotropic AS4/3501-6\n");

    const ff::OrthotropicPly ply = ff::as4_3501_6(125e-6);

    // Layup: [0/45/-45/90/90/-45/45/0] — 8 plies (symmetric)
    const std::vector<double> angles = {0, 45, -45, 90, 90, -45, 45, 0};
    std::vector<ff::LaminatePly> stack;
    for (double a : angles) stack.push_back({ply, a});

    ff::CLTCalculator clt;
    const ff::LaminateABD abd = clt.compute(stack);

    // A₁₁ = A₂₂ for quasi-isotropic laminate
    const double rel_diff = std::abs(abd.A(0,0) - abd.A(1,1)) / abd.A(0,0);
    char detail[256];
    std::snprintf(detail, sizeof(detail),
        "(A11=%.4g, A22=%.4g, rel_diff=%.2e)", abd.A(0,0), abd.A(1,1), rel_diff);
    check(rel_diff < 1e-10, "A11 == A22 (isotropy)", detail);

    // A₁₆ ≈ 0 and A₂₆ ≈ 0  (tolerance: 1 N/m)
    check_abs(abd.A(0,2), 0.0, 1.0, "A16 ≈ 0 [N/m]");
    check_abs(abd.A(1,2), 0.0, 1.0, "A26 ≈ 0 [N/m]");

    // B = 0 for symmetric layup (tolerance: 1 mN)
    check_abs(abd.B.norm(), 0.0, 1e-3, "B == 0 [N]");
}

// ── Case 3: ISA 1976 standard atmosphere ──────────────────────────────────────
// Reference: docs/test_cases.md §Case 3 — NOAA-S/T 76-1562 (1976).
// Tolerance: ±0.1% on all quantities.

static void test_isa_atmosphere() {
    std::printf("\nCase 3 — ISA 1976 atmosphere\n");

    // h = 0 m: T = 288.15 K, p = 101325 Pa, ρ = 1.2250 kg/m³, a = 340.29 m/s
    {
        const auto s = ff::isa1976(0.0);
        check_rel(s.temperature,    288.15,  0.001, "T(0m) [K]");
        check_rel(s.pressure,    101325.0,   0.001, "p(0m) [Pa]");
        check_rel(s.density,       1.2250,   0.001, "rho(0m) [kg/m³]");
        check_rel(s.speed_of_sound, 340.29,  0.001, "a(0m) [m/s]");
    }

    // h = 11 000 m: T = 216.65 K, ρ = 0.3639 kg/m³, a = 295.15 m/s
    {
        const auto s = ff::isa1976(11000.0);
        check_rel(s.temperature,   216.65,  0.001, "T(11km) [K]");
        check_rel(s.density,       0.3639,  0.001, "rho(11km) [kg/m³]");
        check_rel(s.speed_of_sound, 295.15, 0.001, "a(11km) [m/s]");
    }

    // h = 20 000 m: T = 216.65 K (isothermal), ρ = 0.0880 kg/m³
    {
        const auto s = ff::isa1976(20000.0);
        check_rel(s.temperature, 216.65, 0.001, "T(20km) [K]");
        check_rel(s.density,     0.0880, 0.001, "rho(20km) [kg/m³]");
    }
}

// ── Case 4: VLM rectangular flat plate, AR = 5, α = 5° ──────────────────────
// Reference: docs/test_cases.md §Case 4 — Prandtl lifting-line theory.
// Expected: CL ∈ [0.36, 0.42]  (±10% of analytical CL = 0.392).
// AIC matrix: finite and non-singular.

static void test_vlm_rectangular() {
    std::printf("\nCase 4 — VLM rectangular flat plate AR=5, alpha=5deg\n");

    // Rectangular fin: span = 0.5 m, chord = 0.1 m → AR = 2·0.5²/(0.5·0.1) = 10
    // To get AR = 5: use semi-span so that 2·b²/S = 5
    // → b = 0.5 m, c = 0.2 m  → S = 0.5 × 0.2 = 0.1 m², AR = 2×0.25/0.1 = 5
    ff::FinGeometry geo;
    geo.span         = 0.5;
    geo.root_chord   = 0.2;
    geo.tip_chord    = 0.2;  // rectangular: no taper
    geo.sweep_length = 0.0;
    geo.thickness    = 0.001;

    const ff::FlightCondition cond = ff::FlightCondition::from_velocity(0.0, 50.0);

    ff::VortexLattice vlm;
    const double alpha_rad = 5.0 * M_PI / 180.0;
    ff::VLMOptions opts{8, 12, alpha_rad};
    const auto result = vlm.solve(geo, cond, opts);

    // CL should be positive and within ±10% of Prandtl analytical value 0.392
    check(result.CL > 0.0, "CL > 0");
    check(result.CL >= 0.36 && result.CL <= 0.42, "CL in [0.36, 0.42]",
          (std::string("CL=") + std::to_string(result.CL)).c_str());

    // AIC matrix must be non-singular (solution must converge)
    const int n = static_cast<int>(result.panels.size());
    check(n == 8 * 12, "panel count == 96");
    check(result.aic.rows() == n && result.aic.cols() == n, "AIC is square");

    // Check that Γ vector is finite (no NaN/Inf)
    bool gamma_finite = true;
    for (double g : result.gamma)
        if (!std::isfinite(g)) { gamma_finite = false; break; }
    check(gamma_finite, "all gamma values are finite");
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    std::printf("=== fin_flutter C++ benchmark validation ===\n");

    test_clt_unidirectional();
    test_clt_quasi_isotropic();
    test_isa_atmosphere();
    test_vlm_rectangular();

    std::printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}
