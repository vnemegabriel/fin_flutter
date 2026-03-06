import '../../core/materials/orthotropic_ply.dart';
import '../../core/materials/clt_calculator.dart';
import '../../core/materials/laminate.dart';
import '../../core/models/fin_geometry.dart';

/// Bounds for a design variable.
class DesignBound {
  final double min;
  final double max;
  const DesignBound(this.min, this.max);

  double normalize(double value) =>
      ((value - min) / (max - min)).clamp(0.0, 1.0);
  double denormalize(double normalized) =>
      min + normalized.clamp(0.0, 1.0) * (max - min);
}

/// Design variable bounds for optimization.
class DesignBounds {
  final DesignBound span;
  final DesignBound rootChord;
  final DesignBound tipChord;
  final DesignBound sweepAngle; // degrees
  final DesignBound plyAngle;   // degrees, per ply
  final DesignBound thickness;  // per ply, meters

  const DesignBounds({
    this.span = const DesignBound(0.05, 0.5),
    this.rootChord = const DesignBound(0.05, 0.4),
    this.tipChord = const DesignBound(0.0, 0.35),
    this.sweepAngle = const DesignBound(0.0, 60.0),
    this.plyAngle = const DesignBound(-90.0, 90.0),
    this.thickness = const DesignBound(50e-6, 500e-6),
  });
}

/// Design variable vector for optimization.
class DesignVariables {
  double span;         // m
  double rootChord;    // m
  double tipChord;     // m
  double sweepAngle;   // degrees

  List<double> plyAngles;       // degrees per ply
  List<double> plyThicknesses;  // m per ply
  String materialName;          // material database key

  DesignVariables({
    required this.span,
    required this.rootChord,
    required this.tipChord,
    required this.sweepAngle,
    required this.plyAngles,
    required this.plyThicknesses,
    required this.materialName,
  });

  int get plyCount => plyAngles.length;

  double get totalThickness => plyThicknesses.fold(0.0, (s, t) => s + t);

  /// Encode to normalized [0,1] vector for optimizer.
  List<double> toNormalized(DesignBounds bounds) {
    final v = <double>[
      bounds.span.normalize(span),
      bounds.rootChord.normalize(rootChord),
      bounds.tipChord.normalize(tipChord),
      bounds.sweepAngle.normalize(sweepAngle),
    ];
    for (var i = 0; i < plyAngles.length; i++) {
      v.add(bounds.plyAngle.normalize(plyAngles[i]));
      v.add(bounds.thickness.normalize(plyThicknesses[i]));
    }
    return v;
  }

  /// Decode from normalized vector.
  factory DesignVariables.fromNormalized(
      List<double> x, DesignBounds bounds, int nPlies, String materialName) {
    final span = bounds.span.denormalize(x[0]);
    final rootChord = bounds.rootChord.denormalize(x[1]);
    final tipChord = bounds.tipChord
        .denormalize(x[2])
        .clamp(0.0, rootChord); // tipChord <= rootChord
    final sweepAngle = bounds.sweepAngle.denormalize(x[3]);

    final angles = <double>[];
    final thicknesses = <double>[];
    for (var i = 0; i < nPlies; i++) {
      final base = 4 + i * 2;
      if (base + 1 < x.length) {
        angles.add(bounds.plyAngle.denormalize(x[base]));
        thicknesses.add(bounds.thickness.denormalize(x[base + 1]));
      } else {
        angles.add(0.0);
        thicknesses.add(125e-6);
      }
    }

    return DesignVariables(
      span: span,
      rootChord: rootChord,
      tipChord: tipChord,
      sweepAngle: sweepAngle,
      plyAngles: angles,
      plyThicknesses: thicknesses,
      materialName: materialName,
    );
  }

  /// Convert to FinGeometry.
  FinGeometry toFinGeometry() {
    final sweepRad = sweepAngle * 3.14159265358979 / 180.0;
    return FinGeometry(
      span: span,
      rootChord: rootChord,
      tipChord: tipChord,
      sweepLength: span * 0.4 * (sweepAngle / 60.0), // Approximate sweep length
      thickness: totalThickness,
    );
  }

  /// Convert to Laminate.
  Laminate toLaminate(OrthotropicPly ply) {
    final plies = <LaminatePly>[];
    for (var i = 0; i < plyCount; i++) {
      final p = OrthotropicPly(
        name: ply.name,
        E1: ply.E1,
        E2: ply.E2,
        G12: ply.G12,
        nu12: ply.nu12,
        rho: ply.rho,
        t: plyThicknesses[i],
      );
      plies.add(LaminatePly(p, plyAngles[i]));
    }
    return Laminate(plies: plies);
  }

  @override
  String toString() =>
      'DesignVariables(span=${span.toStringAsFixed(3)}m, '
      'rc=${rootChord.toStringAsFixed(3)}m, '
      'plies=$plyCount)';
}
