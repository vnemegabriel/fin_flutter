import '../../math/matrix.dart';
import '../../math/sparse_matrix.dart';
import '../../materials/clt_calculator.dart';
import '../element/kirchhoff_element.dart';
import '../element/mindlin_element.dart';
import '../mesh/boundary_condition.dart';
import '../mesh/mesh.dart';

/// Assembles global stiffness [K] and mass [M] matrices from element contributions.
class GlobalAssembler {
  /// Assemble global matrices for Kirchhoff elements.
  ///
  /// [mesh] - FEA mesh
  /// [abd] - laminate ABD matrices from CLT
  /// [bc] - boundary conditions (clamped DOFs)
  /// Returns (K_global, M_global) as dense matrices.
  ({Matrix K, Matrix M}) assembleKirchhoff({
    required Mesh mesh,
    required LaminateABD abd,
    required BoundaryConditions bc,
  }) {
    final n = mesh.totalDof;
    final K = Matrix(n, n);
    final M = Matrix(n, n);
    final D = abd.D;  // Flexural stiffness matrix [3×3]
    final rhot = abd.areaWeight; // kg/m²

    final element = KirchhoffElement();

    for (final elem in mesh.elements) {
      final coords = _elementCoords(elem, mesh);
      final Ke = element.computeStiffness(coords, D);
      final Me = element.computeMass(coords, rhot);
      final dofs = KirchhoffElement.globalDofIndices(elem.nodeIds);

      _assembleElement(K, Ke, dofs);
      _assembleElement(M, Me, dofs);
    }

    _applyBoundaryConditions(K, M, bc, n);
    return (K: K, M: M);
  }

  /// Assemble global matrices for Mindlin elements.
  ({Matrix K, Matrix M}) assembleMindlin({
    required Mesh mesh,
    required LaminateABD abd,
    required BoundaryConditions bc,
  }) {
    final n = mesh.totalDof;
    final K = Matrix(n, n);
    final M = Matrix(n, n);
    final D = abd.D;
    final rhot = abd.areaWeight;

    // Shear stiffness: Ds = kappa * [A44, 0; 0, A55]
    // Approximate using isotropic shear: A44 = A55 = kappa * h * G
    // From CLT: A66 = G12 * h for 0-degree laminate
    // Use simplified: Ds_diag = kappa * A.get(2,2) / ???
    // For now use 5/6 * A.get(2,2) (shear correction)
    const kappa = 5.0 / 6.0;
    final G_eff = abd.A.get(2, 2) / abd.totalThickness; // effective shear modulus
    final Gkappa = kappa * G_eff * abd.totalThickness;
    final Ds = Matrix.fromList([
      [Gkappa, 0.0],
      [0.0, Gkappa],
    ]);

    final element = MindlinMITC4Element();

    for (final elem in mesh.elements) {
      final coords = _elementCoords(elem, mesh);
      final Ke = element.computeStiffness(coords, D, Ds);
      final Me = element.computeMass(coords, rhot);
      final dofs = MindlinMITC4Element.globalDofIndices(elem.nodeIds);

      _assembleElement(K, Ke, dofs);
      _assembleElement(M, Me, dofs);
    }

    _applyBoundaryConditions(K, M, bc, n);
    return (K: K, M: M);
  }

  /// Get physical coordinates of element corner nodes.
  List<List<double>> _elementCoords(QuadElement elem, Mesh mesh) {
    return elem.nodeIds.map((id) {
      final node = mesh.nodeById(id);
      return [node.x, node.y];
    }).toList();
  }

  /// Scatter element matrix into global matrix.
  void _assembleElement(Matrix K, Matrix Ke, List<int> dofs) {
    final n = dofs.length;
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        K.set(dofs[i], dofs[j], K.get(dofs[i], dofs[j]) + Ke.get(i, j));
      }
    }
  }

  /// Apply homogeneous Dirichlet BC: penalize fixed DOFs in K,
  /// and zero out fixed rows/columns in M.
  void _applyBoundaryConditions(
      Matrix K, Matrix M, BoundaryConditions bc, int n) {
    const bigNumber = 1e30;
    for (final dof in bc.fixedDofs) {
      if (dof >= n) continue;
      // Zero the row and column, put large number on diagonal
      for (var j = 0; j < n; j++) {
        K.set(dof, j, 0.0);
        K.set(j, dof, 0.0);
        M.set(dof, j, 0.0);
        M.set(j, dof, 0.0);
      }
      K.set(dof, dof, bigNumber);
      M.set(dof, dof, 1.0); // Non-zero to avoid singularity in generalized eigenproblem
    }
  }
}
