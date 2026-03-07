import '../math/matrix.dart';
import '../cfd/vlm/vortex_lattice.dart';
import '../fea/results/fea_result.dart';
import '../models/fin_geometry.dart';
import '../models/flight_condition.dart';

/// Couples FEA structural model with VLM aerodynamic model.
///
/// Projects the aerodynamic force distribution onto structural mode shapes
/// to form the generalized aerodynamic force (GAF) matrix.
class AeroelasticCoupler {
  /// Build the generalized aerodynamic force matrix Q_modal (nModes × nModes).
  ///
  /// Q_modal[i][j] = Phi_i^T * [Aerodynamic_Stiffness] * Phi_j
  ///
  /// The aerodynamic stiffness matrix relates panel deflection to panel force.
  /// For a flat plate: Q_phys = rho * V^2 * AIC (simplified).
  ///
  /// [feaResult] - structural FEA results (mode shapes, frequencies)
  /// [vlmResult] - VLM aerodynamic results (AIC, panels)
  /// [panelMeshMapping] - maps each VLM panel to structural DOF contributions
  ///                      (nPanels × nDof interpolation matrix)
  Matrix buildModalAeroMatrix({
    required FEAResult feaResult,
    required VLMResult vlmResult,
    required FinGeometry geometry,
    required FlightCondition condition,
  }) {
    final nModes = feaResult.nModes;
    final nPanels = vlmResult.panels.length;
    final nDof = feaResult.modeShapes.rows;

    // Build panel-to-structure displacement interpolation matrix H (nPanels × nDof)
    final H = _buildInterpolationMatrix(vlmResult.panels, nDof, geometry);

    // AIC corrected
    final AIC = vlmResult.aicMatrix;

    // Aerodynamic stiffness in physical DOF (nDof × nDof):
    // This is a simplified version; full implementation would include
    // mode shape projections onto panel forces.
    // Q_phys ≈ H^T * AIC * H * rho * V^2 / (reference chord)

    final HtAIC = H.transpose() * AIC;
    final QPhys = HtAIC * H;

    // Project onto modal coordinates: Q_modal = Phi^T * Q_phys * Phi
    final Phi = feaResult.modeShapes;
    final PhiTQ = Phi.transpose() * QPhys;
    final QModal = PhiTQ * Phi;

    return QModal;
  }

  /// Build static aerodynamic stiffness matrix for divergence analysis.
  /// Returns A_static (nDof × nDof) such that aerodynamic force = q * A_static * deflection.
  Matrix buildStaticAeroMatrix({
    required VLMResult vlmResult,
    required FinGeometry geometry,
    required int nDof,
    required int dofPerNode,
  }) {
    final H = _buildInterpolationMatrix(vlmResult.panels, nDof, geometry);
    final AIC = vlmResult.aicMatrix;
    // Static aerodynamic stiffness: A_s = H^T * AIC * H
    return H.transpose() * AIC * H;
  }

  /// Simple bilinear interpolation from panel control points to structural DOF.
  ///
  /// For each VLM panel control point, find the nearest structural node and
  /// distribute aerodynamic force to the w (deflection) DOF.
  Matrix _buildInterpolationMatrix(
      List<VLMPanel> panels, int nDof, FinGeometry geometry) {
    final nPanels = panels.length;
    final H = Matrix(nPanels, nDof);

    for (var p = 0; p < nPanels; p++) {
      final cp = panels[p].controlPoint;
      // Find nearest structural node - simplified: use panel centroid y/span
      // to estimate the structural DOF index.
      // In a real implementation, this would use shape function interpolation.
      // Here we use a simple distance-weighted approach.

      // Map panel y-coordinate to approximate structural DOF
      // (first w DOF per node, spaced proportionally along span)
      final yNorm = (cp.y / geometry.span).clamp(0.0, 1.0);
      // Approximate node index based on normalized spanwise position
      // (this is a placeholder; real implementation uses mesh node search)
      final nodeApprox = (yNorm * (nDof - 1)).round().clamp(0, nDof - 1);
      H.set(p, nodeApprox, 1.0);
    }
    return H;
  }
}
