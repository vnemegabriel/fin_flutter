import '../../core/fea/solver/flutter_solver.dart';
import '../../core/fea/solver/divergence_solver.dart';
import '../../core/fea/results/fea_result.dart';
import '../../core/cfd/cfd_engine.dart';
import '../../core/aeroelastic/stability_margin.dart';

/// Complete output from a flutter analysis run.
class FlutterAnalysisOutput {
  final FEAResult feaResult;
  final CFDResult cfdResult;
  final FlutterResult flutterResult;
  final DivergenceResult divergenceResult;
  final double maxFlightVelocity;

  const FlutterAnalysisOutput({
    required this.feaResult,
    required this.cfdResult,
    required this.flutterResult,
    required this.divergenceResult,
    required this.maxFlightVelocity,
  });

  /// Flutter margin (V_flutter / V_max_flight).
  double? get flutterMargin => flutterResult.flutterSpeed != null
      ? StabilityMargin.flutterMargin(flutterResult.flutterSpeed!, maxFlightVelocity)
      : null;

  /// Divergence margin.
  double? get divergenceMargin => divergenceResult.vDivergence != null
      ? StabilityMargin.divergenceMargin(divergenceResult.vDivergence!, maxFlightVelocity)
      : null;

  /// Flutter safety status.
  SafetyStatus? get flutterStatus =>
      flutterMargin != null ? StabilityMargin.flutterStatus(flutterMargin!) : null;

  @override
  String toString() =>
      'FlutterAnalysisOutput('
      'V_F=${flutterResult.flutterSpeed?.toStringAsFixed(1) ?? "N/A"}m/s, '
      'V_D=${divergenceResult.vDivergence?.toStringAsFixed(1) ?? "N/A"}m/s'
      ')';
}
