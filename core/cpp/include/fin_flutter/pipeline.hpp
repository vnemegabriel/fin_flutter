#pragma once
/// @file pipeline.hpp
/// @brief Top-level analysis pipeline: assembles all modules end-to-end.
///
/// Pipeline steps:
///   1. Structural eigenvalue solve → mode shapes Φ, natural frequencies
///   2. VLM aerodynamics → AIC matrix, panel forces
///   3. Aeroelastic coupling → Q_modal (GAF matrix)
///   4. Flutter U-g sweep → V_flutter
///   5. Static divergence sweep → V_divergence
///
/// Reference: docs/theory.md §8 — Analysis Pipeline.

#include <vector>
#include <stdexcept>
#include <Eigen/Dense>

#include "math/types.hpp"
#include "models/fin_geometry.hpp"
#include "models/flight_condition.hpp"
#include "materials/clt_calculator.hpp"
#include "cfd/vortex_lattice.hpp"
#include "fea/eigenvalue_solver.hpp"
#include "aeroelastic/aeroelastic_coupler.hpp"
#include "aeroelastic/flutter_solver.hpp"

namespace ff {

/// @brief Input to the full analysis pipeline.
struct PipelineInput {
    FinGeometry    geometry;
    LaminateABD    laminate;    ///< Pre-computed CLT result.
    FlightCondition condition;

    int n_modes           = 6;  ///< Number of structural modes to extract.
    int chordwise_panels  = 4;  ///< VLM chordwise panel count.
    int spanwise_panels   = 8;  ///< VLM spanwise panel count.
    double alpha_rad      = 2.0 * M_PI / 180.0; ///< Angle of attack [rad].

    double v_min = 10.0;   ///< Velocity sweep start [m/s].
    double v_max = 800.0;  ///< Velocity sweep end   [m/s].
    int    v_steps = 200;  ///< Velocity sweep steps.
};

/// @brief Output of the full analysis pipeline.
struct PipelineOutput {
    EigenResult      eigen;       ///< Structural mode shapes and frequencies.
    VLMResult        vlm;         ///< VLM aerodynamic solution.
    FlutterResult    flutter;     ///< V-g diagram and flutter speed.
    DivergenceResult divergence;  ///< Static divergence speed.
};

/// @brief Run the full flutter / divergence analysis pipeline.
///
/// @param input PipelineInput with geometry, laminate, condition, and solver settings.
/// @return PipelineOutput with all results.
inline PipelineOutput run_pipeline(const PipelineInput& inp) {
    PipelineOutput out;

    // ── 1. Structural model (placeholder beam model from CLT D-matrix) ────────
    //
    // For a cantilever fin plate, a simplified beam model gives approximate
    // natural frequencies. Full FEA with Kirchhoff/Mindlin elements is
    // provided by the FEAEngine in the Dart reference implementation.
    // Here we build a diagonal K/M from Euler-Bernoulli beam theory as
    // a stand-in until the C++ FEA elements are implemented.
    //
    // Eq. 4.1 — docs/theory.md: ω_i = β_i² · √(EI / (ρA L⁴))
    // where β_i L = {1.875, 4.694, 7.855, …} for cantilever beam.
    {
        const int nm = inp.n_modes;
        const double D11  = inp.laminate.D11();    // Nm
        const double rho_s = inp.laminate.area_weight; // kg/m²
        const double b     = inp.geometry.span;
        const double c_r   = inp.geometry.root_chord;
        const double S     = inp.geometry.planform_area();

        // Cantilever beam eigenvalues: β_i·L for modes 1…n
        // Source: Timoshenko & Young (1955) Table 5-3
        const std::vector<double> beta_L = {1.8751, 4.6941, 7.8548, 10.9955,
                                            14.1372, 17.2788};

        MatrixXd K = MatrixXd::Zero(nm, nm);
        MatrixXd M = MatrixXd::Zero(nm, nm);

        for (int i = 0; i < nm; ++i) {
            const double bL = (i < static_cast<int>(beta_L.size()))
                              ? beta_L[i] : M_PI * (i + 0.5);
            // Generalised stiffness and mass (M-normalised → M_ii = 1, K_ii = ω²)
            const double omega_i = (bL * bL / (b * b)) * std::sqrt(D11 / (rho_s + 1e-14));
            K(i,i) = omega_i * omega_i;
            M(i,i) = 1.0;
        }

        EigenvalueSolver solver;
        out.eigen = solver.solve(K, M, nm);
    }

    // ── 2. VLM aerodynamics ───────────────────────────────────────────────────
    VortexLattice vlm;
    out.vlm = vlm.solve(inp.geometry, inp.condition,
                        inp.chordwise_panels, inp.spanwise_panels,
                        inp.alpha_rad);

    // ── 3. Aeroelastic coupling ───────────────────────────────────────────────
    AeroelasticCoupler coupler;
    const MatrixXd Q_modal = coupler.build_modal_aero_matrix(
        out.eigen, out.vlm, inp.geometry);

    // ── 4 & 5. Build modal K and M; sweep for flutter and divergence ──────────
    const int nm = out.eigen.n_modes;
    MatrixXd K_modal = MatrixXd::Zero(nm, nm);
    MatrixXd M_modal = MatrixXd::Zero(nm, nm);
    for (int i = 0; i < nm; ++i) {
        K_modal(i,i) = out.eigen.eigenvalues(i); // ω²  (M-normalised → M=I)
        M_modal(i,i) = 1.0;
    }

    // Velocity sweep
    std::vector<double> vels(inp.v_steps);
    const double dv = (inp.v_max - inp.v_min) / (inp.v_steps - 1);
    for (int i = 0; i < inp.v_steps; ++i)
        vels[i] = inp.v_min + i * dv;

    FlutterSolver fl_solver;
    out.flutter = fl_solver.solve_ug(K_modal, M_modal, Q_modal,
                                     inp.condition.density, vels);

    const int n_dof = static_cast<int>(out.eigen.eigenvectors.rows());
    const MatrixXd A_static = coupler.build_static_aero_matrix(
        out.vlm, inp.geometry, n_dof);
    const MatrixXd& Phi = out.eigen.eigenvectors;
    const MatrixXd A_static_modal = Phi.transpose() * A_static * Phi;

    DivergenceSolver div_solver;
    out.divergence = div_solver.solve(K_modal, A_static_modal,
                                      inp.condition.density,
                                      inp.v_min, inp.v_max);

    return out;
}

} // namespace ff
