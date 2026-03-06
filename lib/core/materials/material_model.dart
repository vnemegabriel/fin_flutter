import '../math/matrix.dart';

/// Abstract base for material models.
abstract class MaterialModel {
  String get name;
  double get density; // kg/m³

  /// Returns the 3×3 reduced stiffness matrix Q (or Q̄ transformed).
  Matrix stiffnessMatrix();
}
