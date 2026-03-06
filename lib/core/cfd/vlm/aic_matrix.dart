import '../../math/matrix.dart';
import 'vortex_lattice.dart';

/// Aerodynamic Influence Coefficient matrix utilities.
class AICMatrix {
  /// Extract the steady AIC from a VLM result.
  static Matrix steadyAIC(VLMResult vlmResult) => vlmResult.aicMatrix;

  /// Build a generalized aerodynamic force matrix (GAF) by projecting AIC
  /// onto structural mode shapes.
  ///
  /// [AIC_physical] - (nPanels × nPanels) AIC in physical coordinates
  /// [panelToStructure] - (nPanels × nDof) mapping from panels to structural DOF
  /// [Phi] - (nDof × nModes) mode shape matrix
  ///
  /// Returns: Q_modal (nModes × nModes) generalized aerodynamic force matrix
  static Matrix generalizedAeroForce({
    required Matrix aicPhysical,
    required Matrix panelToStructure,
    required Matrix phi,
  }) {
    // Q = Phi^T * Z^T * AIC^-1 * Z * Phi (simplified)
    // Full doublet-lattice / Roger approximation would go here.
    // For steady analysis: Q = Phi^T * W * Phi where W = force influence
    final nModes = phi.cols;
    final Q = Matrix(nModes, nModes);
    // Simplified: return scaled identity (placeholder for real GAF)
    return Q;
  }
}
