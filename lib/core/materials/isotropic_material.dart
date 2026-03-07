import '../math/matrix.dart';
import 'material_model.dart';

/// Isotropic material (metals, etc.).
class IsotropicMaterial implements MaterialModel {
  @override
  final String name;
  final double E;    // Young's modulus (Pa)
  final double nu;   // Poisson's ratio
  @override
  final double density; // kg/m³

  const IsotropicMaterial({
    required this.name,
    required this.E,
    required this.nu,
    required this.density,
  });

  double get G => E / (2.0 * (1.0 + nu));

  @override
  Matrix stiffnessMatrix() {
    final Q11 = E / (1.0 - nu * nu);
    final Q12 = nu * E / (1.0 - nu * nu);
    final Q66 = G;
    return Matrix.fromList([
      [Q11, Q12, 0.0],
      [Q12, Q11, 0.0],
      [0.0, 0.0, Q66],
    ]);
  }

  // Common aerospace materials
  static const aluminum2024T3 = IsotropicMaterial(
    name: 'Aluminum 2024-T3',
    E: 73.1e9,
    nu: 0.33,
    density: 2780.0,
  );

  static const steel4340 = IsotropicMaterial(
    name: 'Steel 4340',
    E: 200.0e9,
    nu: 0.29,
    density: 7850.0,
  );

  static const titaniumTi6Al4V = IsotropicMaterial(
    name: 'Titanium Ti-6Al-4V',
    E: 114.0e9,
    nu: 0.34,
    density: 4430.0,
  );
}
