#pragma once
/// @file flutter_solver.hpp
/// @brief U-g (k-method) aeroelastic flutter solver and divergence solver.
///
/// Reference: docs/theory.md §7 — Aeroelastic Stability Analysis.
/// Bisplinghoff, Ashley & Halfman (1955) "Aeroelasticity", §5.4–5.5.
///
/// U-g formulation:
///   The system matrix at velocity V and dynamic pressure q = ½ρV² is:
///     Z = M⁻¹ · (K − q · Q_modal)
///   For each mode i, the damping indicator g_i is:
///     g_i ≈ q · Q_ii / K_ii − 1       (quasi-steady, diagonal approximation)
///   Flutter occurs when g_i crosses zero from below.

#include <vector>
#include <cmath>
#include <optional>
#include <Eigen/Dense>
#include "../math/types.hpp"

namespace ff {

/// @brief One point on the V-g or V-f diagram.
struct VGPoint {
    double velocity;  ///< True airspeed V [m/s].
    double frequency; ///< Modal frequency f [Hz].
    double damping;   ///< Structural damping indicator g [-] (flutter when g ≥ 0).
    int    mode;      ///< Mode index.
};

/// @brief Flutter analysis output.
struct FlutterResult {
    std::optional<double> flutter_speed;     ///< Flutter speed V_F [m/s], if found.
    std::optional<int>    flutter_mode;      ///< Mode that first goes unstable.
    std::optional<double> flutter_frequency; ///< Frequency at V_F [Hz].

    std::vector<std::vector<VGPoint>> vg_curves; ///< V-g curves per mode.
    std::vector<std::vector<VGPoint>> vf_curves; ///< V-f curves per mode.

    bool has_flutter()    const { return flutter_speed.has_value(); }
};

/// @brief Divergence analysis output.
struct DivergenceResult {
    std::optional<double> divergence_speed; ///< Static divergence speed V_D [m/s].
};

/// @brief U-g flutter solver: sweeps velocity and locates the g=0 crossing.
class FlutterSolver {
public:
    /// @brief Perform U-g velocity sweep.
    ///
    /// Algorithm — Eq. 7.3 — Bisplinghoff et al. (1955):
    ///   For each V in sweep:
    ///     q = ½ρV²
    ///     g_i = q · Q_ii / K_ii − 1   (quasi-steady diagonal approximation)
    ///     Flutter detected: g_i changes sign from − to +.
    ///
    /// @param K_modal   Modal stiffness (diagonal, n × n) [N/m per mode].
    /// @param M_modal   Modal mass (diagonal, n × n) [kg per mode].
    /// @param Q_modal   Generalized aerodynamic force matrix (n × n).
    /// @param rho       Air density [kg/m³].
    /// @param velocities Velocity sweep [m/s].
    /// @return FlutterResult with V-g curves and flutter speed.
    FlutterResult solve_ug(const MatrixXd& K_modal,
                           const MatrixXd& M_modal,
                           const MatrixXd& Q_modal,
                           double rho,
                           const std::vector<double>& velocities) const
    {
        const int n_modes = static_cast<int>(K_modal.rows());
        const int n_vel   = static_cast<int>(velocities.size());

        FlutterResult res;
        res.vg_curves.resize(n_modes);
        res.vf_curves.resize(n_modes);

        std::vector<double> prev_g(n_modes, -1.0); // initialise stable

        for (int vi = 0; vi < n_vel; ++vi) {
            const double V = velocities[vi];
            const double q = 0.5 * rho * V * V;

            for (int m = 0; m < n_modes; ++m) {
                const double K_ii = K_modal(m, m);
                const double M_ii = M_modal(m, m);
                const double Q_ii = Q_modal(m, m);

                // Natural frequency at this airspeed (quasi-steady shift)
                double omega2 = (M_ii > 1e-14)
                    ? (K_ii - q * Q_ii) / M_ii : 0.0;
                double freq = (omega2 > 0) ? std::sqrt(omega2) / (2.0 * M_PI) : 0.0;

                // Damping indicator: g < 0 → stable, g ≥ 0 → flutter
                // Eq. 7.3 — Bisplinghoff et al. (1955) §5.4
                double g = (K_ii > 1e-14) ? q * Q_ii / K_ii - 1.0 : 0.0;

                res.vg_curves[m].push_back({V, freq, g, m});
                res.vf_curves[m].push_back({V, freq, g, m});

                // Detect g crossing zero from below
                if (!res.flutter_speed && vi > 0
                    && prev_g[m] < 0.0 && g >= 0.0) {
                    const double V_prev = velocities[vi - 1];
                    const double g_prev = prev_g[m];
                    // Linear interpolation — Eq. 7.4
                    const double V_F = V_prev + (V - V_prev) * (-g_prev) / (g - g_prev);
                    res.flutter_speed     = V_F;
                    res.flutter_mode      = m;
                    res.flutter_frequency = freq;
                }
                prev_g[m] = g;
            }
        }
        return res;
    }
};

/// @brief Static aeroelastic divergence solver.
///
/// Divergence speed V_D is the lowest velocity for which
/// the aeroelastic stiffness matrix is singular:
///   det(K_modal − q · A_static_modal) = 0
///
/// Reference: docs/theory.md §7.2 — Bisplinghoff et al. (1955) §8.
class DivergenceSolver {
public:
    /// @brief Find static divergence speed by sweeping velocity and detecting determinant sign change.
    ///
    /// Optimization (Phase 4.2): Use determinant sign crossing instead of full eigensolve
    /// for ~5-10× speedup. Divergence occurs when det(K_ae) = 0.
    ///
    /// @param K_modal        Modal stiffness (n × n).
    /// @param A_static_modal Projected static aero stiffness (n × n).
    /// @param rho            Air density [kg/m³].
    /// @param v_min          Minimum velocity [m/s].
    /// @param v_max          Maximum velocity [m/s].
    /// @param n_steps        Number of velocity steps.
    /// @return DivergenceResult with divergence speed if found.
    DivergenceResult solve(const MatrixXd& K_modal,
                           const MatrixXd& A_static_modal,
                           double rho,
                           double v_min = 1.0,
                           double v_max = 1500.0,
                           int    n_steps = 500) const
    {
        DivergenceResult res;
        const double dv = (v_max - v_min) / (n_steps - 1);

        double prev_det = 1.0; // Will track sign of determinant

        for (int i = 0; i < n_steps; ++i) {
            const double V = v_min + i * dv;
            const double q = 0.5 * rho * V * V;

            // Aeroelastic stiffness matrix K_ae = K − q·A_s
            const MatrixXd K_ae = K_modal - q * A_static_modal;

            // Divergence: det(K_ae) crosses zero (cheaper than full eigensolve)
            // det = 0 when matrix becomes singular (lowest eigenvalue = 0)
            const double det = K_ae.determinant();

            // Detect sign change: divergence found when determinant crosses zero
            if (i > 0 && prev_det > 0.0 && det <= 0.0) {
                res.divergence_speed = V;
                break;
            }
            prev_det = det;
        }
        return res;
    }
};

} // namespace ff
