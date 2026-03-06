import 'dart:typed_data';
import '../../core/models/fin_geometry.dart';
import '../../core/models/flight_condition.dart';
import 'ork_geometry.dart';
import 'ork_parser.dart';

/// High-level facade for importing OpenRocket .ork files.
class OrkImporter {
  final OrkParser _parser = OrkParser();

  /// Import from raw file bytes.
  Future<OrkImportResult> import(Uint8List fileBytes) async {
    final doc = await _parser.parse(fileBytes);

    final finGeo = doc.primaryFinSet?.toFinGeometry();
    final flightProfile = doc.primaryFlightProfile;

    // Determine worst-case flight condition from flight profile
    FlightCondition? worstCase;
    if (flightProfile != null && flightProfile.dataPoints.isNotEmpty) {
      // Find point of maximum dynamic pressure
      var maxQPoint = flightProfile.dataPoints.first;
      for (final pt in flightProfile.dataPoints) {
        final rho = FlightCondition.isa(altitude: pt.altitude, mach: 0).density;
        final q = 0.5 * rho * pt.velocity * pt.velocity;
        final qMax = 0.5 *
            FlightCondition.isa(altitude: maxQPoint.altitude, mach: 0).density *
            maxQPoint.velocity *
            maxQPoint.velocity;
        if (q > qMax) maxQPoint = pt;
      }
      worstCase = FlightCondition.fromVelocity(
        altitude: maxQPoint.altitude,
        velocity: maxQPoint.velocity,
      );
    }

    return OrkImportResult(
      document: doc,
      finGeometry: finGeo,
      flightCondition: worstCase,
      flightProfile: flightProfile,
    );
  }
}

/// Result of importing an OpenRocket file.
class OrkImportResult {
  final OrkDocument document;
  final FinGeometry? finGeometry;
  final FlightCondition? flightCondition;
  final OrkFlightProfile? flightProfile;

  const OrkImportResult({
    required this.document,
    this.finGeometry,
    this.flightCondition,
    this.flightProfile,
  });

  bool get hasGeometry => finGeometry != null;
  bool get hasFlightCondition => flightCondition != null;
  String get rocketName => document.rocketName;
}
