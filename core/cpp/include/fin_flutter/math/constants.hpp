#pragma once
/// @file constants.hpp
/// @brief Physical constants, mathematical constants, and numerical tolerances.
///
/// Consolidates constants scattered throughout the codebase for centralized
/// maintenance and documentation.

namespace ff {

// ──────────────────────────────────────────────────────────────────────────────
// ISA 1976 Standard Atmosphere Constants
// Reference: NACA Report 1235, docs/theory.md §1
// ──────────────────────────────────────────────────────────────────────────────

/// @brief Sea-level temperature T₀ [K].
static constexpr double kISA_T0 = 288.15;

/// @brief Sea-level pressure p₀ [Pa].
static constexpr double kISA_P0 = 101325.0;

/// @brief Specific gas constant for dry air R [J/(kg·K)].
static constexpr double kISA_R = 287.058;

/// @brief Gravitational acceleration g₀ [m/s²].
static constexpr double kISA_G0 = 9.80665;

/// @brief Heat capacity ratio γ = Cp/Cv [-] (diatomic air).
static constexpr double kISA_GAMMA = 1.4;

/// @brief Tropospheric lapse rate L [K/m].
static constexpr double kISA_L_TROP = 0.0065;

/// @brief Tropopause transition altitude [m].
static constexpr double kISA_H_TROP = 11000.0;

/// @brief Stratosphere 1 upper limit [m].
static constexpr double kISA_H_STRAT1 = 20000.0;

/// @brief Stratosphere 2 upper limit [m].
static constexpr double kISA_H_STRAT2 = 32000.0;

// ──────────────────────────────────────────────────────────────────────────────
// Mathematical Constants
// ──────────────────────────────────────────────────────────────────────────────

/// @brief π (pi).
static constexpr double kPI = 3.14159265358979323846;

/// @brief 2π (two pi).
static constexpr double k2PI = 2.0 * kPI;

/// @brief 1/(4π) — common factor in Biot-Savart law.
static constexpr double kBiot_Savart_Coeff = 1.0 / (4.0 * kPI);

// ──────────────────────────────────────────────────────────────────────────────
// Aerodynamic & VLM Constants
// ──────────────────────────────────────────────────────────────────────────────

/// @brief Cutoff radius to avoid singularity at vortex core [m].
/// Used in Biot-Savart finite and semi-infinite segment computations.
/// Reference: Katz & Plotkin (2001) §10.
static constexpr double kVortex_Cutoff = 1e-6;

/// @brief Vortex strength sign convention:
/// - Bound vortex (A→B): Γ > 0 (counterclockwise viewed from +z)
/// - Right trailing (B→+∞): Γ > 0 (downwash)
/// - Left trailing (A→+∞): Γ > 0 (downwash)
/// This produces AIC(i,i) < 0 (downwash diagonal) → CL > 0 for α > 0.
/// Reference: docs/architecture.md §5.3 (horseshoe vortex sign convention)

// ──────────────────────────────────────────────────────────────────────────────
// Classical Lamination Theory (CLT) Constants
// ──────────────────────────────────────────────────────────────────────────────

/// @brief Laminate midplane location convention: z = 0 at midplane.
/// Plies span z ∈ [−h/2, +h/2] where h is total thickness.
/// Reference: Reddy (2004) §3.3, docs/theory.md §2.2.

// ──────────────────────────────────────────────────────────────────────────────
// Test Tolerances & Acceptance Criteria
// ──────────────────────────────────────────────────────────────────────────────

namespace test_tolerances {

/// @brief CLT matrix relative tolerance (±0.2%).
/// Applied to A₁₁, A₁₂, A₆₆, D₁₁ per Reddy (2004).
static constexpr double kCLT_Rel = 0.002;

/// @brief CLT matrix absolute tolerance for near-zero elements [N/m].
/// Applied to B matrix (symmetric laminates should have B ≈ 0).
static constexpr double kCLT_Abs = 1e-3;

/// @brief ISA 1976 relative tolerance (±0.5%).
/// Applied to temperature, pressure, density, speed of sound.
static constexpr double kISA_Rel = 0.005;

/// @brief VLM lift coefficient (CL) empirical acceptable range [lower bound].
/// Reference: docs/test_cases.md Case 4 (under investigation for magnitude).
static constexpr double kVLM_CL_Min = 0.36;

/// @brief VLM lift coefficient (CL) empirical acceptable range [upper bound].
static constexpr double kVLM_CL_Max = 0.42;

} // namespace test_tolerances

} // namespace ff
