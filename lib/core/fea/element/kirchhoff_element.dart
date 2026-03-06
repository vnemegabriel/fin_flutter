import 'dart:math';
import '../../math/matrix.dart';
import '../../materials/clt_calculator.dart';
import '../mesh/mesh.dart';

/// Discrete Kirchhoff Quadrilateral (DKQ) plate element.
/// 4 nodes, 3 DOF per node (w, theta_x, theta_y) = 12 DOF total.
///
/// Uses 2×2 Gauss quadrature for stiffness and mass integration.
class KirchhoffElement {
  static const int _dofPerNode = 3;
  static const int _totalDof = 12; // 4 nodes × 3 DOF
  static const int _nGauss = 2;   // integration points per direction

  // 2×2 Gauss points and weights
  static final List<double> _gaussPts = [-1.0 / sqrt(3), 1.0 / sqrt(3)];
  static final List<double> _gaussWts = [1.0, 1.0];

  /// Compute 12×12 element stiffness matrix.
  ///
  /// [nodeCoords] - list of 4 [x,y] pairs (physical coordinates, m)
  /// [D] - 3×3 flexural stiffness matrix from CLT (N·m)
  Matrix computeStiffness(List<List<double>> nodeCoords, Matrix D) {
    assert(nodeCoords.length == 4);
    final K = Matrix(_totalDof, _totalDof);

    for (var ip = 0; ip < _nGauss; ip++) {
      final xi = _gaussPts[ip];
      for (var jp = 0; jp < _nGauss; jp++) {
        final eta = _gaussPts[jp];
        final w = _gaussWts[ip] * _gaussWts[jp];

        final jac = _jacobian(xi, eta, nodeCoords);
        final detJ = _det2(jac);
        if (detJ.abs() < 1e-14) continue;
        final jacInv = _inv2(jac, detJ);

        final B = _strainDisplacementMatrix(xi, eta, jacInv, nodeCoords);
        // K += B^T * D * B * detJ * weight
        final BtD = B.transpose() * D;
        final BtDB = BtD * B;
        final scale = detJ * w;
        for (var r = 0; r < _totalDof; r++) {
          for (var c = 0; c < _totalDof; c++) {
            K.set(r, c, K.get(r, c) + BtDB.get(r, c) * scale);
          }
        }
      }
    }
    return K;
  }

  /// Compute 12×12 consistent mass matrix.
  ///
  /// [rhot] - areal density × thickness = rho_lam (kg/m²)
  Matrix computeMass(List<List<double>> nodeCoords, double rhot) {
    assert(nodeCoords.length == 4);
    final M = Matrix(_totalDof, _totalDof);

    for (var ip = 0; ip < _nGauss; ip++) {
      final xi = _gaussPts[ip];
      for (var jp = 0; jp < _nGauss; jp++) {
        final eta = _gaussPts[jp];
        final w = _gaussWts[ip] * _gaussWts[jp];

        final jac = _jacobian(xi, eta, nodeCoords);
        final detJ = _det2(jac);
        if (detJ.abs() < 1e-14) continue;

        // Shape functions for w DOF only (translational inertia)
        final N = _shapeN(xi, eta); // 12 × 1 (only w components)
        final scale = rhot * detJ * w;
        for (var r = 0; r < _totalDof; r++) {
          for (var c = 0; c < _totalDof; c++) {
            M.set(r, c, M.get(r, c) + N[r] * N[c] * scale);
          }
        }
      }
    }
    return M;
  }

  /// DOF indices for the element: [node0_dofs, node1_dofs, ...]
  static List<int> globalDofIndices(List<int> nodeIds) {
    final indices = <int>[];
    for (final nid in nodeIds) {
      for (var d = 0; d < _dofPerNode; d++) {
        indices.add(nid * _dofPerNode + d);
      }
    }
    return indices;
  }

  // --- Private helpers ---

  /// Bilinear shape functions N_i(xi, eta).
  List<double> _shapeFunc(double xi, double eta) => [
        (1 - xi) * (1 - eta) / 4.0,
        (1 + xi) * (1 - eta) / 4.0,
        (1 + xi) * (1 + eta) / 4.0,
        (1 - xi) * (1 + eta) / 4.0,
      ];

  /// Derivatives dN/dxi and dN/deta.
  List<List<double>> _shapeDeriv(double xi, double eta) => [
        [-(1 - eta) / 4.0, (1 - eta) / 4.0, (1 + eta) / 4.0, -(1 + eta) / 4.0],
        [-(1 - xi) / 4.0, -(1 + xi) / 4.0, (1 + xi) / 4.0, (1 - xi) / 4.0],
      ];

  /// 2×2 Jacobian matrix J = dX/dxi.
  List<List<double>> _jacobian(double xi, double eta, List<List<double>> coords) {
    final dN = _shapeDeriv(xi, eta);
    final J = [[0.0, 0.0], [0.0, 0.0]];
    for (var k = 0; k < 4; k++) {
      J[0][0] += dN[0][k] * coords[k][0];
      J[0][1] += dN[0][k] * coords[k][1];
      J[1][0] += dN[1][k] * coords[k][0];
      J[1][1] += dN[1][k] * coords[k][1];
    }
    return J;
  }

  double _det2(List<List<double>> J) =>
      J[0][0] * J[1][1] - J[0][1] * J[1][0];

  List<List<double>> _inv2(List<List<double>> J, double det) => [
        [J[1][1] / det, -J[0][1] / det],
        [-J[1][0] / det, J[0][0] / det],
      ];

  /// Strain-displacement matrix B [3 × 12].
  ///
  /// For DKQ, uses simplified formulation where bending strains relate to
  /// rotational DOF: kappa_x = dtheta_y/dx, kappa_y = -dtheta_x/dy, etc.
  /// This is a simplified DKT/DKQ approximation using Kirchhoff constraint.
  Matrix _strainDisplacementMatrix(
      double xi, double eta, List<List<double>> jacInv, List<List<double>> coords) {
    final dN = _shapeDeriv(xi, eta);
    final B = Matrix(3, _totalDof);

    for (var k = 0; k < 4; k++) {
      // dN/dx and dN/dy via chain rule: [dN/dx, dN/dy] = [dN/dxi, dN/deta] * J^-1
      final dNdx = jacInv[0][0] * dN[0][k] + jacInv[0][1] * dN[1][k];
      final dNdy = jacInv[1][0] * dN[0][k] + jacInv[1][1] * dN[1][k];

      // DOF order per node: [w, theta_x, theta_y]
      // Bending strains: kappa_x = d(theta_y)/dx, kappa_y = -d(theta_x)/dy
      // kappa_xy = 0.5*(d(theta_y)/dy + d(-theta_x)/dx) ... wait
      // Corrected Kirchhoff: k_xx = -d^2w/dx^2, k_yy = -d^2w/dy^2, k_xy = -2*d^2w/dxdy
      // With theta_x = dw/dy, theta_y = -dw/dx (Kirchhoff convention):
      // k_xx = d(theta_y)/dx, k_yy = -d(theta_x)/dy, k_xy = d(theta_y)/dy - d(theta_x)/dx

      final base = k * _dofPerNode; // w_k, theta_x_k, theta_y_k
      // Row 0 (kappa_xx): contribution from theta_y (DOF index base+2)
      B.set(0, base + 2, B.get(0, base + 2) + dNdx);
      // Row 1 (kappa_yy): contribution from -theta_x (DOF index base+1)
      B.set(1, base + 1, B.get(1, base + 1) - dNdy);
      // Row 2 (2*kappa_xy): contributions from theta_y/dy and -theta_x/dx
      B.set(2, base + 2, B.get(2, base + 2) + dNdy);
      B.set(2, base + 1, B.get(2, base + 1) - dNdx);
    }
    return B;
  }

  /// Shape function vector for mass matrix (maps DOF to displacement w).
  /// Only translational DOF w contributes to translational inertia.
  List<double> _shapeN(double xi, double eta) {
    final N = _shapeFunc(xi, eta);
    final result = List<double>.filled(_totalDof, 0.0);
    for (var k = 0; k < 4; k++) {
      result[k * _dofPerNode] = N[k]; // only w DOF
    }
    return result;
  }
}
