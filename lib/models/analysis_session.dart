import '../core/models/fin_geometry.dart';
import '../core/models/flight_condition.dart';
import '../core/materials/clt_calculator.dart';
import '../core/materials/laminate.dart';
import '../modules/flutter_analysis/flutter_output.dart';

/// A complete analysis session: geometry + materials + results.
class AnalysisSession {
  final String id;
  final String name;
  final FinGeometry geometry;
  final Laminate laminate;
  final FlightCondition flightCondition;
  final FlutterAnalysisOutput? results;
  final DateTime createdAt;
  final DateTime updatedAt;

  const AnalysisSession({
    required this.id,
    required this.name,
    required this.geometry,
    required this.laminate,
    required this.flightCondition,
    this.results,
    required this.createdAt,
    required this.updatedAt,
  });

  AnalysisSession copyWith({
    String? name,
    FinGeometry? geometry,
    Laminate? laminate,
    FlightCondition? flightCondition,
    FlutterAnalysisOutput? results,
  }) =>
      AnalysisSession(
        id: id,
        name: name ?? this.name,
        geometry: geometry ?? this.geometry,
        laminate: laminate ?? this.laminate,
        flightCondition: flightCondition ?? this.flightCondition,
        results: results ?? this.results,
        createdAt: createdAt,
        updatedAt: DateTime.now(),
      );

  bool get hasResults => results != null;

  @override
  String toString() => 'AnalysisSession($id: $name)';
}
