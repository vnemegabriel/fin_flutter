#pragma once
/// @file orthotropic_ply.hpp
/// @brief Single unidirectional composite ply: elastic constants and stiffness matrices.
///
/// Reference: docs/theory.md §3 — Classical Lamination Theory.
/// Reddy (2004) "Mechanics of Laminated Composite Plates and Shells", §1.3.

#include <cmath>
#include <stdexcept>
#include <string>
#include "../math/types.hpp"

namespace ff {

/// @brief Material properties of a single unidirectional orthotropic ply.
struct OrthotropicPly {
    std::string name;

    // ── Elastic constants (all in Pa) ─────────────────────────────────────────
    double E1;   ///< Longitudinal (fibre-direction) modulus E₁ [Pa].
    double E2;   ///< Transverse modulus E₂ [Pa].
    double G12;  ///< In-plane shear modulus G₁₂ [Pa].
    double nu12; ///< Major Poisson's ratio ν₁₂ [-].
    double rho;  ///< Density ρ [kg/m³].
    double t;    ///< Ply thickness [m].

    // ── Strength properties (Pa, with typical defaults) ───────────────────────
    double Xt  = 1500e6; ///< Longitudinal tensile strength [Pa].
    double Xc  = 1200e6; ///< Longitudinal compressive strength [Pa].
    double Yt  =   50e6; ///< Transverse tensile strength [Pa].
    double Yc  =  200e6; ///< Transverse compressive strength [Pa].
    double S12 =   70e6; ///< In-plane shear strength [Pa].

    // ── Derived constants ─────────────────────────────────────────────────────

    /// @brief Minor Poisson's ratio ν₂₁ = ν₁₂ · E₂/E₁  (reciprocal relation).
    /// Eq. 3.2 — docs/theory.md
    double nu21() const { return nu12 * E2 / E1; }

    /// @brief Reduced stiffness matrix [Q] in principal material axes (3×3).
    ///
    /// Voigt notation: [Q₁₁ Q₁₂ 0; Q₁₂ Q₂₂ 0; 0 0 Q₆₆].
    ///
    /// Eq. 3.3 — Reddy (2004) Eq. 1.3.71
    ///   Q₁₁ = E₁ / (1 − ν₁₂ν₂₁)
    ///   Q₂₂ = E₂ / (1 − ν₁₂ν₂₁)
    ///   Q₁₂ = ν₁₂·E₂ / (1 − ν₁₂ν₂₁)
    ///   Q₆₆ = G₁₂
    ///
    /// @throws std::invalid_argument if material properties are outside physical bounds.
    Matrix3d reduced_stiffness() const {
        // Input validation: ensure physical plausibility
        if (E1 <= 0.0 || E2 <= 0.0 || G12 <= 0.0) {
            throw std::invalid_argument(
                "OrthotropicPly::reduced_stiffness: moduli must be positive (E1, E2, G12 > 0)");
        }
        if (nu12 < 0.0 || nu12 >= 0.5) {
            throw std::invalid_argument(
                "OrthotropicPly::reduced_stiffness: Poisson's ratio must be in [0, 0.5)");
        }

        const double delta = 1.0 - nu12 * nu21();
        const double Q11 = E1  / delta;
        const double Q22 = E2  / delta;
        const double Q12 = nu12 * E2 / delta;
        const double Q66 = G12;
        Matrix3d Q;
        Q << Q11, Q12, 0.0,
             Q12, Q22, 0.0,
             0.0, 0.0, Q66;
        return Q;
    }

    /// @brief Transformed reduced stiffness [Q̄] at ply angle θ [degrees].
    ///
    /// Uses the standard tensor rotation — Eq. 3.4 — Reddy (2004) Eq. 1.3.81.
    ///   Q̄₁₁ = Q₁₁m⁴ + 2(Q₁₂+2Q₆₆)m²n² + Q₂₂n⁴
    ///   Q̄₁₂ = (Q₁₁+Q₂₂−4Q₆₆)m²n² + Q₁₂(m⁴+n⁴)
    ///   (etc.)
    ///
    /// @param theta_deg Ply orientation angle θ [degrees] (0 = fibre along x-axis).
    Matrix3d transformed_stiffness(double theta_deg) const {
        const double theta = theta_deg * M_PI / 180.0;
        const double m = std::cos(theta);
        const double n = std::sin(theta);
        const double m2 = m * m,  n2 = n * n;
        const double mn = m * n,  m2n2 = m2 * n2;

        const Matrix3d Q = reduced_stiffness();
        const double Q11 = Q(0,0), Q22 = Q(1,1), Q12 = Q(0,1), Q66 = Q(2,2);

        const double Qb11 = Q11*m2*m2 + 2.0*(Q12+2.0*Q66)*m2n2 + Q22*n2*n2;
        const double Qb12 = (Q11+Q22-4.0*Q66)*m2n2 + Q12*(m2*m2+n2*n2);
        const double Qb22 = Q11*n2*n2 + 2.0*(Q12+2.0*Q66)*m2n2 + Q22*m2*m2;
        const double Qb16 = (Q11-Q12-2.0*Q66)*m2*mn - (Q22-Q12-2.0*Q66)*mn*n2;
        const double Qb26 = (Q11-Q12-2.0*Q66)*mn*n2 - (Q22-Q12-2.0*Q66)*m2*mn;
        const double Qb66 = (Q11+Q22-2.0*Q12-2.0*Q66)*m2n2 + Q66*(m2*m2+n2*n2);

        Matrix3d Qb;
        Qb << Qb11, Qb12, Qb16,
              Qb12, Qb22, Qb26,
              Qb16, Qb26, Qb66;
        return Qb;
    }
};

// ── Common material presets ───────────────────────────────────────────────────

/// @brief AS4/3501-6 carbon/epoxy ply (t = 125 µm).
/// Source: Reddy (2004) Table 3.3; Tsai & Hahn (1980).
inline OrthotropicPly as4_3501_6(double t_ply = 125e-6) {
    return {"AS4/3501-6",
            142e9, 10.3e9, 7.2e9, 0.27,
            1580.0, t_ply};
}

/// @brief T300/5208 carbon/epoxy ply (t = 125 µm).
/// Source: Tsai & Hahn (1980).
inline OrthotropicPly t300_5208(double t_ply = 125e-6) {
    return {"T300/5208",
            181e9, 10.3e9, 7.17e9, 0.28,
            1600.0, t_ply};
}

/// @brief IM7/8552 carbon/epoxy ply (t = 131 µm).
inline OrthotropicPly im7_8552(double t_ply = 131e-6) {
    return {"IM7/8552",
            171e9, 9.08e9, 5.29e9, 0.32,
            1570.0, t_ply};
}

/// @brief E-glass/Epoxy ply (t = 200 µm).
inline OrthotropicPly eglass_epoxy(double t_ply = 200e-6) {
    return {"E-glass/Epoxy",
            45.6e9, 16.2e9, 5.83e9, 0.278,
            2100.0, t_ply};
}

} // namespace ff
