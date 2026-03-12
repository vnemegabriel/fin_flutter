#pragma once
/// @file fin_geometry.hpp
/// @brief Trapezoidal fin planform geometry.
///
/// Reference: docs/theory.md §2 — Fin Geometry Definitions.

#include <cmath>
#include <stdexcept>

namespace ff {

/// @brief Trapezoidal fin geometry definition.
///
/// Coordinate system: x along chord (root LE at origin), y along span.
/// All dimensions in SI (metres).
struct FinGeometry {
    double span;         ///< Semi-span b (m) — root to tip.
    double root_chord;   ///< Root chord c_r (m).
    double tip_chord;    ///< Tip chord  c_t (m).
    double sweep_length; ///< Leading-edge sweep offset x_s (m).
    double thickness;    ///< Fin thickness t (m).
    int    fin_count{4}; ///< Number of fins on the rocket.

    // ── Derived quantities ────────────────────────────────────────────────────

    /// @brief Planform area S = ½(c_r + c_t)·b  [m²].
    /// Eq. 2.1 — docs/theory.md §2.2
    double planform_area() const {
        return 0.5 * (root_chord + tip_chord) * span;
    }

    /// @brief Aspect ratio AR = b² / S (using full-span b = 2·semi-span).
    /// Eq. 2.2 — docs/theory.md §2.3
    double aspect_ratio() const {
        const double S = planform_area();
        return (S > 0) ? (2.0 * span * span) / S : 0.0;
    }

    /// @brief Leading-edge sweep angle Λ [rad].
    /// Eq. 2.3 — docs/theory.md §2.4
    double sweep_angle_rad() const { return std::atan2(sweep_length, span); }

    /// @brief Leading-edge sweep angle [deg].
    double sweep_angle_deg() const { return sweep_angle_rad() * 180.0 / M_PI; }

    /// @brief Taper ratio λ = c_t / c_r.
    double taper_ratio() const {
        return (root_chord > 0) ? tip_chord / root_chord : 0.0;
    }

    /// @brief Mean aerodynamic chord MAC [m].
    /// Eq. 2.5 — docs/theory.md §2.5
    double mac() const {
        const double lambda = taper_ratio();
        return (2.0/3.0) * root_chord * (1.0 + lambda + lambda*lambda) / (1.0 + lambda);
    }

    /// @brief Leading-edge x-coordinate at spanwise station y [m].
    double leading_edge_x(double y) const {
        return sweep_length * (y / span);
    }

    /// @brief Chord length at spanwise station y [m].
    double chord_at(double y) const {
        return root_chord + (tip_chord - root_chord) * (y / span);
    }

    /// @brief Trailing-edge x-coordinate at spanwise station y [m].
    double trailing_edge_x(double y) const {
        return leading_edge_x(y) + chord_at(y);
    }
};

} // namespace ff
