import 'clt_calculator.dart';
import 'orthotropic_ply.dart';

/// A composite laminate defined by an ordered ply stack.
class Laminate {
  final List<LaminatePly> plies;
  final String name;

  const Laminate({required this.plies, this.name = 'Custom Laminate'});

  /// Total laminate thickness (m).
  double get totalThickness => plies.fold(0.0, (sum, p) => sum + p.ply.t);

  /// Total areal density (kg/m²).
  double get arealDensity => plies.fold(0.0, (sum, p) => sum + p.ply.rho * p.ply.t);

  /// Number of plies.
  int get plyCount => plies.length;

  /// Compute ABD matrices.
  LaminateABD computeABD() => CLTCalculator().compute(plies);

  /// Check if laminate is symmetric about midplane.
  bool get isSymmetric {
    final n = plies.length;
    for (var i = 0; i < n ~/ 2; i++) {
      final top = plies[i];
      final bot = plies[n - 1 - i];
      if (top.ply.name != bot.ply.name ||
          top.angleDegrees != bot.angleDegrees ||
          top.ply.t != bot.ply.t) {
        return false;
      }
    }
    return true;
  }

  /// Check if laminate is balanced (every +θ ply has a -θ ply).
  bool get isBalanced {
    final angles = <double, int>{};
    for (final p in plies) {
      angles[p.angleDegrees] = (angles[p.angleDegrees] ?? 0) + 1;
    }
    for (final entry in angles.entries) {
      final angle = entry.key;
      if (angle == 0.0 || angle == 90.0 || angle == -90.0) continue;
      final count = entry.value;
      final complementCount = angles[-angle] ?? 0;
      if (count != complementCount) return false;
    }
    return true;
  }

  /// Standard laminates.
  factory Laminate.unidirectional(OrthotropicPly ply, int nPlies) {
    return Laminate(
      name: '[0]_${nPlies}',
      plies: List.generate(nPlies, (_) => LaminatePly(ply, 0.0)),
    );
  }

  factory Laminate.crossPly(OrthotropicPly ply, int nPairsEach) {
    final stack = <LaminatePly>[];
    for (var i = 0; i < nPairsEach; i++) {
      stack.add(LaminatePly(ply, 0.0));
      stack.add(LaminatePly(ply, 90.0));
    }
    return Laminate(name: '[0/90]_${nPairsEach}', plies: stack);
  }

  factory Laminate.quasiIsotropic(OrthotropicPly ply) {
    return Laminate(
      name: '[0/+45/-45/90]_s',
      plies: [
        LaminatePly(ply, 0.0),
        LaminatePly(ply, 45.0),
        LaminatePly(ply, -45.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, -45.0),
        LaminatePly(ply, 45.0),
        LaminatePly(ply, 0.0),
      ],
    );
  }

  @override
  String toString() => 'Laminate($name, ${plyCount} plies, h=${(totalThickness*1000).toStringAsFixed(2)}mm)';
}
