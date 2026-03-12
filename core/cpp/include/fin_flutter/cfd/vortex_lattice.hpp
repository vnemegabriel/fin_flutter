#pragma once
/// @file vortex_lattice.hpp
/// @brief Vortex Lattice Method (VLM) solver for fin aerodynamics.
///
/// Reference: docs/theory.md §5 — Vortex Lattice Method.
/// Katz & Plotkin (2001) "Low-Speed Aerodynamics", §12.
///
/// Discretisation: each panel carries a horseshoe vortex.
///   Bound vortex: at 1/4-chord line.
///   Control point: at 3/4-chord (Kutta condition via Pistolesi theorem).
///   Normal: z-direction (flat-plate assumption).
///   Wake: semi-infinite, trailing in +x direction.

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "biot_savart.hpp"
#include "../math/vec3.hpp"
#include "../models/fin_geometry.hpp"
#include "../models/flight_condition.hpp"

namespace ff {

/// @brief A single VLM panel (horseshoe vortex element).
struct VLMPanel {
    int    id;
    Vec3   bound_a;      ///< Inboard end of bound vortex (1/4-chord).
    Vec3   bound_b;      ///< Outboard end of bound vortex.
    Vec3   control_pt;   ///< Control point (3/4-chord).
    Vec3   normal;       ///< Surface normal at control point.
    double area;         ///< Panel area [m²].
    double span_width;   ///< Spanwise panel width Δy [m].
};

/// @brief Output of the VLM solve.
struct VLMResult {
    std::vector<double>     gamma;     ///< Vortex strength Γ per panel [m²/s].
    std::vector<double>     panel_cl;  ///< Lift coefficient per panel.
    double                  CL;        ///< Total lift coefficient.
    double                  CDi;       ///< Induced drag coefficient.
    Eigen::MatrixXd         aic;       ///< AIC matrix (nPanels × nPanels).
    std::vector<VLMPanel>   panels;    ///< Panel geometry.
};

/// @brief Vortex Lattice Method solver.
class VortexLattice {
public:
    /// @brief Solve VLM for a trapezoidal fin.
    ///
    /// @param geometry           Fin planform geometry.
    /// @param condition          Flight condition (V∞, ρ).
    /// @param chordwise_panels   Number of panels in chord direction n_x.
    /// @param spanwise_panels    Number of panels in span  direction n_y.
    /// @param alpha_rad          Angle of attack [rad].
    /// @return VLMResult with Γ distribution, lift, drag, AIC.
    VLMResult solve(const FinGeometry&     geometry,
                    const FlightCondition& condition,
                    int    chordwise_panels,
                    int    spanwise_panels,
                    double alpha_rad) const
    {
        const int nx = chordwise_panels;
        const int ny = spanwise_panels;
        const int n  = nx * ny;

        const auto panels = build_panels(geometry, nx, ny);
        const auto AIC    = build_aic(panels, condition.velocity);

        // RHS: − (V∞ · n̂) for each panel — Eq. 5.8 (flow-tangency BC)
        Eigen::VectorXd rhs(n);
        const Vec3 v_inf{condition.velocity * std::cos(alpha_rad),
                         0.0,
                         condition.velocity * std::sin(alpha_rad)};
        for (int i = 0; i < n; ++i)
            rhs(i) = -panels[i].normal.dot(v_inf);

        // Solve AIC · Γ = rhs  (dense LU)
        Eigen::VectorXd gamma_vec = AIC.lu().solve(rhs);

        // ── Forces via Kutta-Joukowski — Eq. 5.9 ─────────────────────────────
        const double q   = condition.dynamic_pressure();
        const double S   = geometry.planform_area();
        std::vector<double> gamma(n), panel_cl(n);
        double total_lift = 0.0;

        for (int i = 0; i < n; ++i) {
            gamma[i] = gamma_vec(i);
            const double L  = condition.density * condition.velocity
                              * gamma[i] * panels[i].span_width;
            panel_cl[i] = (q * panels[i].area > 0) ? L / (q * panels[i].area) : 0.0;
            total_lift += L;
        }

        const double CL = (q * S > 0) ? total_lift / (q * S) : 0.0;

        // Induced drag (simplified Trefftz-plane summation)
        double CDi = 0.0;
        for (int i = 0; i < n; ++i) {
            const double w_i = induced_normal_velocity(panels, gamma, panels[i].control_pt);
            CDi += gamma[i] * w_i;
        }
        CDi = std::abs(CDi) / (2.0 * condition.velocity * S);

        return {gamma, panel_cl, CL, CDi, AIC, panels};
    }

    // ── Accessors (for testing and coupling) ─────────────────────────────────

    /// @brief Build panel geometry for a fin.
    std::vector<VLMPanel> build_panels(const FinGeometry& geo, int nx, int ny) const {
        std::vector<VLMPanel> panels;
        panels.reserve(nx * ny);
        int id = 0;

        for (int j = 0; j < ny; ++j) {
            const double y0   = static_cast<double>(j)   / ny * geo.span;
            const double y1   = static_cast<double>(j+1) / ny * geo.span;
            const double y_mid = (y0 + y1) / 2.0;

            const double x_le_mid   = geo.leading_edge_x(y_mid);
            const double chord_mid  = geo.chord_at(y_mid);

            for (int i = 0; i < nx; ++i) {
                const double xi0 = static_cast<double>(i)   / nx;
                const double xi1 = static_cast<double>(i+1) / nx;

                // Bound vortex at 1/4-chord within the panel chord extent
                const double x_bound = x_le_mid + (xi0 + (xi1-xi0)/4.0) * chord_mid;
                // Control point at 3/4-chord (Pistolesi theorem)
                const double x_ctrl  = x_le_mid + (xi0 + (xi1-xi0)*3.0/4.0) * chord_mid;

                panels.push_back({
                    id++,
                    Vec3{x_bound, y0, 0.0},
                    Vec3{x_bound, y1, 0.0},
                    Vec3{x_ctrl,  y_mid, 0.0},
                    Vec3{0.0, 0.0, 1.0}, // flat-plate normal
                    chord_mid * (xi1-xi0) * (y1-y0),
                    y1 - y0
                });
            }
        }
        return panels;
    }

    /// @brief Build Aerodynamic Influence Coefficient matrix.
    ///
    /// AIC[i,j] = normal velocity at control point i due to unit Γ on panel j.
    /// Eq. 5.7 — Katz & Plotkin (2001) §12.3.
    Eigen::MatrixXd build_aic(const std::vector<VLMPanel>& panels,
                               double /*V_inf*/) const {
        const int n = static_cast<int>(panels.size());
        Eigen::MatrixXd AIC(n, n);
        const Vec3 wake{1.0, 0.0, 0.0}; // wake trailing downstream (+x)

        for (int i = 0; i < n; ++i) {
            const Vec3& cp = panels[i].control_pt;
            const Vec3& ni = panels[i].normal;

            for (int j = 0; j < n; ++j) {
                const Vec3& A = panels[j].bound_a;
                const Vec3& B = panels[j].bound_b;

                const Vec3 v_bound  = BiotSavart::finite_segment(cp, A, B, 1.0);
                const Vec3 v_trail_a = BiotSavart::semi_infinite(cp, A,  1.0, wake);  // Same sign as bound/B trailing
                const Vec3 v_trail_b = BiotSavart::semi_infinite(cp, B,  1.0, wake);
                const Vec3 v_total  = v_bound + v_trail_a + v_trail_b;

                AIC(i, j) = ni.dot(v_total);
            }
        }
        return AIC;
    }

private:
    /// @brief z-component of total induced velocity at point P from all panels.
    double induced_normal_velocity(const std::vector<VLMPanel>& panels,
                                   const std::vector<double>& gamma,
                                   const Vec3& P) const {
        const Vec3 wake{1.0, 0.0, 0.0};
        Vec3 v{};
        for (int j = 0; j < static_cast<int>(panels.size()); ++j) {
            const Vec3& A = panels[j].bound_a;
            const Vec3& B = panels[j].bound_b;
            const double g = gamma[j];
            v = v + BiotSavart::finite_segment(P, A, B,  g)
                  + BiotSavart::semi_infinite(P, A,  g, wake)  // Same sign as B trailing
                  + BiotSavart::semi_infinite(P, B,  g, wake);
        }
        return v.z;
    }
};

} // namespace ff
