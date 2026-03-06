import '../../core/models/fin_geometry.dart';
import '../../core/models/flight_condition.dart';
import '../../core/materials/clt_calculator.dart';
import '../../core/fea/fea_engine.dart';
import '../../core/cfd/cfd_engine.dart';

/// Input for a full flutter/divergence analysis.
class FlutterAnalysisInput {
  final FinGeometry geometry;
  final LaminateABD laminate;
  final FlightCondition flightCondition;
  final FEAConfig feaConfig;
  final CFDConfig cfdConfig;
  final double velocityMin;
  final double velocityMax;
  final int velocitySteps;
  final double maxFlightVelocity;

  const FlutterAnalysisInput({
    required this.geometry,
    required this.laminate,
    required this.flightCondition,
    this.feaConfig = const FEAConfig(),
    this.cfdConfig = const CFDConfig(),
    this.velocityMin = 10.0,
    this.velocityMax = 600.0,
    this.velocitySteps = 50,
    this.maxFlightVelocity = 300.0,
  });

  List<double> get velocitySweep {
    final step = (velocityMax - velocityMin) / (velocitySteps - 1);
    return List.generate(velocitySteps, (i) => velocityMin + i * step);
  }
}
