#pragma once
/// @file aeroelastic_coupler.hpp
/// @brief Projects VLM aerodynamic forces onto structural mode shapes.
///
/// Reference: docs/theory.md §6 — Aeroelastic Coupling.
///
/// The generalized aerodynamic force (GAF) matrix is:
///   Q_modal = Φᵀ · H^T · [AIC] · H · Φ
/// where H (nPanels × nDOF) maps structural deflections to panel displacements.

#include <Eigen/Dense>
#include "../math/types.hpp"
#include "../cfd/vortex_lattice.hpp"
#include "../fea/eigenvalue_solver.hpp"
#include "../models/fin_geometry.hpp"

namespace ff {

/// @brief Couples FEA mode shapes with VLM aerodynamic influence coefficients.
class AeroelasticCoupler {
public:
    /// @brief Build the generalized aerodynamic force matrix Q_modal (nModes × nModes).
    ///
    /// Eq. 6.3 — docs/theory.md:
    ///   Q_modal[i,j] = Φᵢᵀ · Q_phys · Φⱼ
    ///   Q_phys       = H^T · AIC · H
    ///
    /// @param eigen_res  Structural eigenproblem result (mode shapes Φ, n_DOF × n_modes).
    /// @param vlm_result VLM result (AIC, panels).
    /// @param geometry   Fin geometry (for spanwise interpolation).
    /// @return Q_modal (n_modes × n_modes).
    MatrixXd build_modal_aero_matrix(const EigenResult&  eigen_res,
                                     const VLMResult&    vlm_result,
                                     const FinGeometry&  geometry) const
    {
        const int n_dof   = static_cast<int>(eigen_res.eigenvectors.rows());
        const int n_modes = eigen_res.n_modes;
        const int n_pan   = static_cast<int>(vlm_result.panels.size());

        const MatrixXd H   = build_interpolation(vlm_result.panels, n_dof, geometry);
        const MatrixXd Q_phys = H.transpose() * vlm_result.aic * H; // nDOF × nDOF
        const MatrixXd& Phi   = eigen_res.eigenvectors;              // nDOF × nModes

        return Phi.transpose() * Q_phys * Phi; // nModes × nModes
    }

    /// @brief Build static aerodynamic stiffness matrix for divergence analysis.
    ///
    /// A_static (n_DOF × n_DOF) satisfies: aerodynamic force = q · A_static · u.
    /// Eq. 6.5 — docs/theory.md: A_static = H^T · AIC · H.
    MatrixXd build_static_aero_matrix(const VLMResult&   vlm_result,
                                      const FinGeometry& geometry,
                                      int                n_dof) const
    {
        const MatrixXd H = build_interpolation(vlm_result.panels, n_dof, geometry);
        return H.transpose() * vlm_result.aic * H;
    }

private:
    /// @brief Build panel-to-structure displacement interpolation matrix H.
    ///
    /// H (nPanels × nDOF): H[p, nodeApprox] = 1, all other entries 0.
    ///
    /// Maps each VLM panel control point to the nearest structural DOF
    /// using a simple nearest-node lookup along the span.
    /// A full implementation would use isoparametric shape functions.
    MatrixXd build_interpolation(const std::vector<VLMPanel>& panels,
                                 int n_dof,
                                 const FinGeometry& geo) const
    {
        const int n_pan = static_cast<int>(panels.size());
        MatrixXd H = MatrixXd::Zero(n_pan, n_dof);

        for (int p = 0; p < n_pan; ++p) {
            const double y_norm = std::clamp(panels[p].control_pt.y / geo.span,
                                             0.0, 1.0);
            const int node = static_cast<int>(std::round(y_norm * (n_dof - 1)));
            H(p, node) = 1.0;
        }
        return H;
    }
};

} // namespace ff
