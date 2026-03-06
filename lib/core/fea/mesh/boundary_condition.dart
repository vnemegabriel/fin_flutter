import 'mesh.dart';

/// Represents constrained DOFs in the FEA model.
class BoundaryConditions {
  /// Set of global DOF indices that are fixed (zero displacement).
  final Set<int> fixedDofs;

  const BoundaryConditions({required this.fixedDofs});

  /// Create clamped root boundary condition for a fin.
  ///
  /// All nodes at the root (y ≈ 0) are fully clamped:
  /// w = 0, theta_x = 0, theta_y = 0 (all DOF zeroed).
  factory BoundaryConditions.clampedRoot(Mesh mesh) {
    final fixed = <int>{};
    final dofPerNode = mesh.dofPerNode;

    for (final node in mesh.nodes) {
      if (node.y.abs() < 1e-10) {
        // Root node: fix all DOF
        final base = node.id * dofPerNode;
        for (var d = 0; d < dofPerNode; d++) {
          fixed.add(base + d);
        }
      }
    }
    return BoundaryConditions(fixedDofs: fixed);
  }

  bool isFixed(int dof) => fixedDofs.contains(dof);
}
