#pragma once
/// @file biot_savart.hpp
/// @brief Biot-Savart law for vortex-induced velocities.
///
/// Reference: docs/theory.md §5 — Vortex Lattice Method.
/// Katz & Plotkin (2001) "Low-Speed Aerodynamics", §10.1–10.2.

#include <cmath>
#include "../math/vec3.hpp"

namespace ff {

/// @brief Cutoff radius to avoid singularity at the vortex core.
static constexpr double kVortexCutoff = 1e-6;

/// @brief Biot-Savart induced-velocity kernels for VLM horseshoe vortices.
namespace BiotSavart {

/// @brief Velocity induced by a finite vortex segment A→B with circulation Γ.
///
/// Eq. 5.3 — Katz & Plotkin (2001) Eq. 10.11:
///   V = (Γ/4π) · (r₁×r₂) / |r₁×r₂|² · (B−A)·(r̂₁−r̂₂)
///   NOTE: (B−A) is unnormalized; this preserves the |AB| factor needed so that
///   the formula gives V ∝ 1/h (perpendicular distance), not 1/(h·|AB|).
///
/// @param P   Field point where velocity is evaluated.
/// @param A   Vortex segment start.
/// @param B   Vortex segment end.
/// @param gamma Circulation strength Γ [m²/s].
/// @return Induced velocity vector [m/s].
inline Vec3 finite_segment(const Vec3& P, const Vec3& A, const Vec3& B,
                            double gamma) {
    const Vec3 r1 = P - A;
    const Vec3 r2 = P - B;
    const Vec3 cross = r1.cross(r2);
    const double cross_mag2 = cross.dot(cross);

    if (cross_mag2 < kVortexCutoff * kVortexCutoff)
        return {};

    const double r1_mag = r1.norm();
    const double r2_mag = r2.norm();
    if (r1_mag < kVortexCutoff || r2_mag < kVortexCutoff)
        return {};

    // cosTerm = (B−A) · (r̂₁ − r̂₂) — use UNNORMALIZED (B−A) so the |AB|
    // factor is retained; this makes h = |r₁×r₂|/|AB| the true perpendicular
    // distance and the induced velocity ∝ 1/h as expected.
    const Vec3   AB       = B - A;
    const double cos_term = AB.dot(r1.normalized())
                          - AB.dot(r2.normalized());
    const double strength = gamma / (4.0 * M_PI) * cos_term / cross_mag2;
    return cross * strength;
}

/// @brief Velocity induced by a semi-infinite trailing vortex from point A.
///
/// The vortex extends from A to +∞ along `direction` (wake direction).
///
/// Eq. 5.4 — Katz & Plotkin (2001) Eq. 10.15:
///   V = (Γ/4π) · (r×d̂) / |r×d̂|² · (1 + r·d̂/|r|)
///
/// @param P         Field point.
/// @param A         Vortex origin.
/// @param gamma     Circulation Γ [m²/s].
/// @param direction Unit vector of wake direction (typically +x).
/// @return Induced velocity vector [m/s].
inline Vec3 semi_infinite(const Vec3& P, const Vec3& A, double gamma,
                           const Vec3& direction) {
    const Vec3  r       = P - A;
    const Vec3  r_cross_d = r.cross(direction);
    const double cross_mag2 = r_cross_d.dot(r_cross_d);

    if (cross_mag2 < kVortexCutoff * kVortexCutoff)
        return {};

    const double r_dot_d = r.dot(direction);
    const double r_mag   = r.norm();
    const double factor  = gamma / (4.0 * M_PI) * (1.0 + r_dot_d / r_mag) / cross_mag2;
    return r_cross_d * factor;
}

} // namespace BiotSavart
} // namespace ff
