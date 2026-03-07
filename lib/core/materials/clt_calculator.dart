import '../math/matrix.dart';
import 'orthotropic_ply.dart';

/// A ply in the laminate stack: ply + orientation angle.
class LaminatePly {
  final OrthotropicPly ply;
  final double angleDegrees; // ply angle in degrees

  const LaminatePly(this.ply, this.angleDegrees);
}

/// [A], [B], [D] matrices and derived properties from Classical Lamination Theory.
class LaminateABD {
  final Matrix A; // Extensional stiffness (3×3), N/m
  final Matrix B; // Bending-extension coupling (3×3), N
  final Matrix D; // Flexural stiffness (3×3), N·m

  final double totalThickness; // m
  final double areaWeight;     // kg/m²

  const LaminateABD({
    required this.A,
    required this.B,
    required this.D,
    required this.totalThickness,
    required this.areaWeight,
  });

  /// 6×6 combined ABD matrix.
  Matrix get combinedABD {
    final abd = Matrix(6, 6);
    for (var i = 0; i < 3; i++) {
      for (var j = 0; j < 3; j++) {
        abd.set(i, j, A.get(i, j));
        abd.set(i, j + 3, B.get(i, j));
        abd.set(i + 3, j, B.get(i, j));
        abd.set(i + 3, j + 3, D.get(i, j));
      }
    }
    return abd;
  }

  /// Effective in-plane modulus Ex (Pa).
  double get effectiveEx {
    final h = totalThickness;
    if (h <= 0) return 0;
    final denom = A.get(0, 0) * A.get(1, 1) - A.get(0, 1) * A.get(0, 1);
    return denom / (A.get(1, 1) * h);
  }

  /// Effective bending stiffness D11 (N·m).
  double get effectiveBendingStiffness => D.get(0, 0);

  @override
  String toString() =>
      'LaminateABD(h=${(totalThickness * 1000).toStringAsFixed(2)}mm, '
      'w=${areaWeight.toStringAsFixed(2)}kg/m²)';
}

/// Classical Lamination Theory calculator.
/// Computes [A], [B], [D] matrices for a laminate stack.
class CLTCalculator {
  /// Compute ABD matrices for the given ply stack.
  ///
  /// [stack] is ordered from bottom to top (z increasing upward).
  /// Midplane is at z=0; plies span from z = -h/2 to z = +h/2.
  LaminateABD compute(List<LaminatePly> stack) {
    if (stack.isEmpty) {
      return LaminateABD(
        A: Matrix(3, 3),
        B: Matrix(3, 3),
        D: Matrix(3, 3),
        totalThickness: 0,
        areaWeight: 0,
      );
    }

    final totalH = stack.fold(0.0, (sum, lp) => sum + lp.ply.t);
    final A = Matrix(3, 3);
    final B = Matrix(3, 3);
    final D = Matrix(3, 3);

    // z_k coordinates: start at -totalH/2, step through each ply
    var zBottom = -totalH / 2.0;
    for (final lp in stack) {
      final zTop = zBottom + lp.ply.t;
      final Qbar = lp.ply.transformedStiffness(lp.angleDegrees);

      final dz = zTop - zBottom;           // z_k - z_{k-1}
      final dz2 = zTop * zTop - zBottom * zBottom; // z_k^2 - z_{k-1}^2
      final dz3 = zTop * zTop * zTop - zBottom * zBottom * zBottom; // z_k^3 - z_{k-1}^3

      for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) {
          final q = Qbar.get(i, j);
          A.set(i, j, A.get(i, j) + q * dz);
          B.set(i, j, B.get(i, j) + q * dz2 / 2.0);
          D.set(i, j, D.get(i, j) + q * dz3 / 3.0);
        }
      }
      zBottom = zTop;
    }

    final areaWeight = stack.fold(0.0, (sum, lp) => sum + lp.ply.rho * lp.ply.t);

    return LaminateABD(
      A: A,
      B: B,
      D: D,
      totalThickness: totalH,
      areaWeight: areaWeight,
    );
  }

  /// Convenience method: compute ABD for symmetric balanced laminate.
  /// [halfStack] is the bottom half; automatically mirrored.
  LaminateABD computeSymmetric(List<LaminatePly> halfStack) {
    final fullStack = [...halfStack, ...halfStack.reversed.toList()];
    return compute(fullStack);
  }

  /// Compute engineering constants from ABD matrix.
  EngineeringConstants engineeringConstants(LaminateABD abd) {
    final h = abd.totalThickness;
    if (h <= 0) return EngineeringConstants.zero;
    final a = abd.A.inverse();
    final Ex = 1.0 / (a.get(0, 0) * h);
    final Ey = 1.0 / (a.get(1, 1) * h);
    final Gxy = 1.0 / (a.get(2, 2) * h);
    final nuxy = -a.get(0, 1) / a.get(0, 0);
    return EngineeringConstants(Ex: Ex, Ey: Ey, Gxy: Gxy, nuxy: nuxy);
  }
}

/// Effective engineering constants derived from CLT.
class EngineeringConstants {
  final double Ex;   // Pa
  final double Ey;   // Pa
  final double Gxy;  // Pa
  final double nuxy; // dimensionless

  const EngineeringConstants({
    required this.Ex,
    required this.Ey,
    required this.Gxy,
    required this.nuxy,
  });

  static const zero = EngineeringConstants(Ex: 0, Ey: 0, Gxy: 0, nuxy: 0);

  @override
  String toString() =>
      'EngineeringConstants(Ex=${(Ex/1e9).toStringAsFixed(2)}GPa, '
      'Ey=${(Ey/1e9).toStringAsFixed(2)}GPa, '
      'Gxy=${(Gxy/1e9).toStringAsFixed(2)}GPa, '
      'nuxy=${nuxy.toStringAsFixed(3)})';
}
