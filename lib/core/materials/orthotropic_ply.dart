import 'dart:math';
import '../math/matrix.dart';

/// Single unidirectional composite ply properties.
class OrthotropicPly {
  final String name;
  final double E1;    // Longitudinal modulus (Pa)
  final double E2;    // Transverse modulus (Pa)
  final double G12;   // In-plane shear modulus (Pa)
  final double nu12;  // Major Poisson's ratio
  final double rho;   // Density (kg/m³)
  final double t;     // Ply thickness (m)

  // Strength properties (Pa)
  final double xt;    // Longitudinal tensile strength
  final double xc;    // Longitudinal compressive strength
  final double yt;    // Transverse tensile strength
  final double yc;    // Transverse compressive strength
  final double s12;   // In-plane shear strength

  const OrthotropicPly({
    required this.name,
    required this.E1,
    required this.E2,
    required this.G12,
    required this.nu12,
    required this.rho,
    required this.t,
    this.xt = 1500e6,
    this.xc = 1200e6,
    this.yt = 50e6,
    this.yc = 200e6,
    this.s12 = 70e6,
  });

  /// Minor Poisson's ratio (nu21).
  double get nu21 => nu12 * E2 / E1;

  /// Reduced stiffness matrix [Q] in principal material axes (3×3).
  /// [ Q11 Q12  0  ]
  /// [ Q12 Q22  0  ]
  /// [  0   0  Q66 ]
  Matrix get reducedStiffness {
    final denom = 1.0 - nu12 * nu21;
    final Q11 = E1 / denom;
    final Q22 = E2 / denom;
    final Q12 = nu12 * E2 / denom;
    final Q66 = G12;
    return Matrix.fromList([
      [Q11, Q12, 0.0],
      [Q12, Q22, 0.0],
      [0.0, 0.0, Q66],
    ]);
  }

  /// Transformed reduced stiffness [Q̄] at ply angle theta (degrees).
  /// Uses standard transformation matrix.
  Matrix transformedStiffness(double thetaDegrees) {
    final theta = thetaDegrees * pi / 180.0;
    final m = cos(theta);
    final n = sin(theta);
    final m2 = m * m;
    final n2 = n * n;
    final mn = m * n;
    final m2n2 = m2 * n2;

    final denom = 1.0 - nu12 * nu21;
    final Q11 = E1 / denom;
    final Q22 = E2 / denom;
    final Q12 = nu12 * E2 / denom;
    final Q66 = G12;

    final Qb11 = Q11 * m2 * m2 + 2.0 * (Q12 + 2.0 * Q66) * m2n2 + Q22 * n2 * n2;
    final Qb12 = (Q11 + Q22 - 4.0 * Q66) * m2n2 + Q12 * (m2 * m2 + n2 * n2);
    final Qb22 = Q11 * n2 * n2 + 2.0 * (Q12 + 2.0 * Q66) * m2n2 + Q22 * m2 * m2;
    final Qb16 = (Q11 - Q12 - 2.0 * Q66) * m2 * mn - (Q22 - Q12 - 2.0 * Q66) * mn * n2;
    final Qb26 = (Q11 - Q12 - 2.0 * Q66) * mn * n2 - (Q22 - Q12 - 2.0 * Q66) * m2 * mn;
    final Qb66 = (Q11 + Q22 - 2.0 * Q12 - 2.0 * Q66) * m2n2 + Q66 * (m2 * m2 + n2 * n2);

    return Matrix.fromList([
      [Qb11, Qb12, Qb16],
      [Qb12, Qb22, Qb26],
      [Qb16, Qb26, Qb66],
    ]);
  }

  @override
  String toString() => 'OrthotropicPly($name: E1=${(E1/1e9).toStringAsFixed(1)}GPa)';
}
