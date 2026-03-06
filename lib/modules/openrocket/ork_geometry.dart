import 'dart:math';
import '../../core/models/fin_geometry.dart';

/// Fin geometry data extracted from an OpenRocket file.
class OrkFinGeometry {
  final double span;          // m
  final double rootChord;     // m
  final double tipChord;      // m
  final double sweepLength;   // m (leading edge sweep)
  final double thickness;     // m
  final String crossSection;  // 'square', 'rounded', 'airfoil'
  final String materialName;
  final int finCount;
  final double mountDiameter; // rocket body tube diameter (m)

  const OrkFinGeometry({
    required this.span,
    required this.rootChord,
    required this.tipChord,
    required this.sweepLength,
    required this.thickness,
    this.crossSection = 'square',
    this.materialName = 'Unknown',
    this.finCount = 4,
    this.mountDiameter = 0.1,
  });

  /// Convert to FinGeometry for analysis.
  FinGeometry toFinGeometry() => FinGeometry(
        span: span,
        rootChord: rootChord,
        tipChord: tipChord,
        sweepLength: sweepLength,
        thickness: thickness,
        finCount: finCount,
      );

  double get aspectRatio => (2.0 * span * span) / planformArea;
  double get planformArea => 0.5 * (rootChord + tipChord) * span;
  double get sweepAngleDeg => atan2(sweepLength, span) * 180.0 / pi;
  double get taperRatio => rootChord > 0 ? tipChord / rootChord : 0.0;

  @override
  String toString() =>
      'OrkFinGeometry(span=${span.toStringAsFixed(3)}m, '
      'rc=${rootChord.toStringAsFixed(3)}m, '
      'tc=${tipChord.toStringAsFixed(3)}m)';
}

/// Flight data point extracted from OpenRocket simulation.
class OrkFlightDataPoint {
  final double time;      // s
  final double altitude;  // m
  final double velocity;  // m/s
  final double mach;      // dimensionless

  const OrkFlightDataPoint({
    required this.time,
    required this.altitude,
    required this.velocity,
    required this.mach,
  });
}

/// Flight profile from OpenRocket simulation.
class OrkFlightProfile {
  final List<OrkFlightDataPoint> dataPoints;

  const OrkFlightProfile({required this.dataPoints});

  double get maxVelocity => dataPoints.fold(0.0, (m, p) => m > p.velocity ? m : p.velocity);
  double get maxMach => dataPoints.fold(0.0, (m, p) => m > p.mach ? m : p.mach);
  double get maxAltitude => dataPoints.fold(0.0, (m, p) => m > p.altitude ? m : p.altitude);
}

/// Complete OpenRocket document data.
class OrkDocument {
  final String rocketName;
  final List<OrkFinGeometry> finSets;
  final List<OrkFlightProfile> flightProfiles;

  const OrkDocument({
    required this.rocketName,
    required this.finSets,
    required this.flightProfiles,
  });

  bool get hasFinSets => finSets.isNotEmpty;
  bool get hasFlightData => flightProfiles.isNotEmpty;

  OrkFinGeometry? get primaryFinSet => hasFinSets ? finSets.first : null;
  OrkFlightProfile? get primaryFlightProfile =>
      hasFlightData ? flightProfiles.first : null;
}
