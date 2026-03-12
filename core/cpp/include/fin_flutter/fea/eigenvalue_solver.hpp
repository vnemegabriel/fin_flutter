#pragma once
/// @file eigenvalue_solver.hpp
/// @brief Generalized structural eigenvalue solver: [K]{φ} = λ[M]{φ}.
///
/// Reference: docs/theory.md §4 — Finite Element Analysis.
/// Uses Eigen's GeneralizedSelfAdjointEigenSolver (Cholesky-based).
/// For very large systems, a sparse ARPACK/LOBPCG back-end would be substituted.

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "../math/types.hpp"

namespace ff {

/// @brief Result of a structural eigenvalue solve.
struct EigenResult {
    VectorXd eigenvalues;    ///< Natural frequencies squared ω² [rad²/s²], ascending.
    MatrixXd eigenvectors;   ///< Corresponding mode shapes Φ (columns), M-normalised.
    int      n_modes;        ///< Number of modes extracted.

    /// @brief Natural frequencies f_i [Hz].
    std::vector<double> frequencies_hz() const {
        std::vector<double> f(n_modes);
        for (int i = 0; i < n_modes; ++i)
            f[i] = std::sqrt(std::abs(eigenvalues(i))) / (2.0 * M_PI);
        return f;
    }

    /// @brief Natural frequencies ω_i [rad/s].
    std::vector<double> frequencies_rad() const {
        std::vector<double> w(n_modes);
        for (int i = 0; i < n_modes; ++i)
            w[i] = std::sqrt(std::abs(eigenvalues(i)));
        return w;
    }
};

/// @brief Solves the generalized eigenvalue problem [K]{φ} = λ[M]{φ}.
///
/// Assumptions:
///   - [K] is symmetric positive semi-definite (after boundary conditions).
///   - [M] is symmetric positive definite.
///   - n_modes ≤ matrix dimension.
///
/// Method: Eigen::GeneralizedSelfAdjointEigenSolver (Cholesky decomposition).
/// Complexity: O(n³) — suitable for n ≤ ~5 000 DOF.
class EigenvalueSolver {
public:
    /// @brief Compute the first n_modes eigenpairs of [K]{φ} = λ[M]{φ}.
    ///
    /// @param K       Stiffness matrix (n×n, symmetric PSD).
    /// @param M       Mass matrix (n×n, symmetric PD).
    /// @param n_modes Number of lowest modes to return.
    /// @return EigenResult with eigenvalues sorted ascending and M-normalised eigenvectors.
    EigenResult solve(const MatrixXd& K, const MatrixXd& M, int n_modes) const {
        const int n = static_cast<int>(K.rows());
        const int n_req = std::min(n_modes, n);

        // Eq. 4.8 — docs/theory.md: Cholesky-based generalised solve
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(K, M);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("EigenvalueSolver: generalised eigensolve failed");

        // Eigenvalues are already sorted ascending
        const VectorXd& evals = solver.eigenvalues();
        const MatrixXd& evecs = solver.eigenvectors();

        EigenResult res;
        res.n_modes     = n_req;
        res.eigenvalues = evals.head(n_req);
        res.eigenvectors = evecs.leftCols(n_req);

        return res;
    }
};

} // namespace ff
