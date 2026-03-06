import '../../core/materials/material_database.dart';
import '../../core/models/flight_condition.dart';
import '../../core/fea/fea_engine.dart';
import '../../core/cfd/cfd_engine.dart';
import '../flutter_analysis/flutter_analysis_module.dart';
import '../flutter_analysis/flutter_input.dart';
import 'design_variables.dart';

/// Optimization objective and constraint functions.
class FlutterLossFunction {
  final double w1;      // flutter speed weight
  final double w2;      // divergence speed weight
  final double w3;      // mass penalty weight
  final double vRef;    // reference velocity (m/s) - target minimum flutter speed
  final double mRef;    // reference mass (kg/m²) - reference areal density
  final double mu;      // constraint penalty coefficient

  final FlightCondition flightCondition;
  final FEAConfig feaConfig;
  final CFDConfig cfdConfig;

  const FlutterLossFunction({
    this.w1 = 1.0,
    this.w2 = 0.5,
    this.w3 = 0.2,
    required this.vRef,
    this.mRef = 2.0,
    this.mu = 10.0,
    required this.flightCondition,
    this.feaConfig = const FEAConfig(
      spanwiseElements: 4,
      chordwiseElements: 4,
      nModes: 3,
    ),
    this.cfdConfig = const CFDConfig(
      chordwisePanels: 4,
      spanwisePanels: 6,
    ),
  });

  /// Evaluate loss function for given design variables.
  LossEvaluation evaluate(DesignVariables vars) {
    // Get material
    final ply = MaterialDatabase.findCompositePly(vars.materialName) ??
        MaterialDatabase.as43501;

    // Build geometry and laminate
    final geometry = vars.toFinGeometry();
    final laminate = vars.toLaminate(ply);
    final abd = laminate.computeABD();

    // Run analysis
    final module = FlutterAnalysisModule();
    final output = module.analyze(FlutterAnalysisInput(
      geometry: geometry,
      laminate: abd,
      flightCondition: flightCondition,
      feaConfig: feaConfig,
      cfdConfig: cfdConfig,
      velocityMin: 50.0,
      velocityMax: vRef * 2.0,
      velocitySteps: 20,
      maxFlightVelocity: vRef / 1.2,
    ));

    final vFlutter = output.flutterResult.flutterSpeed ?? vRef * 2.0;
    final vDivergence = output.divergenceResult.vDivergence ?? vRef * 2.0;
    final mass = abd.areaWeight * geometry.planformArea;
    final mNorm = mass / mRef.clamp(0.01, double.infinity);

    // Objective function
    final obj = w1 * (vFlutter / vRef - 1.0) * (vFlutter / vRef - 1.0) +
        w2 * (vDivergence / vRef - 1.0) * (vDivergence / vRef - 1.0) +
        w3 * mNorm;

    // Constraint violations
    final constraints = _evaluateConstraints(vars, vFlutter, mass, abd.totalThickness);
    final penalty = constraints.fold(0.0, (sum, g) {
      final v = g > 0 ? g : 0.0;
      return sum + mu * v * v;
    });

    return LossEvaluation(
      objective: obj,
      penalty: penalty,
      total: obj + penalty,
      flutterSpeed: vFlutter,
      divergenceSpeed: vDivergence,
      mass: mass,
      constraints: constraints,
    );
  }

  /// Constraint functions g_i(x) <= 0 (negative = satisfied).
  List<double> _evaluateConstraints(
      DesignVariables vars, double vFlutter, double mass, double totalThickness) {
    return [
      // g1: flutter speed must exceed 1.2 * vRef
      vRef * 1.2 - vFlutter,
      // g2: span within bounds
      vars.span - 0.5,
      0.05 - vars.span,
      // g3: ply angle change rate (manufacturing)
      ...List.generate(vars.plyAngles.length - 1, (i) {
        final delta = (vars.plyAngles[i + 1] - vars.plyAngles[i]).abs();
        return delta - 45.0; // max 45° change per ply
      }),
      // g4: total thickness
      totalThickness - 0.05,
    ];
  }
}

/// Result of a loss function evaluation.
class LossEvaluation {
  final double objective;
  final double penalty;
  final double total;
  final double flutterSpeed;
  final double divergenceSpeed;
  final double mass;
  final List<double> constraints;

  const LossEvaluation({
    required this.objective,
    required this.penalty,
    required this.total,
    required this.flutterSpeed,
    required this.divergenceSpeed,
    required this.mass,
    required this.constraints,
  });

  bool get feasible => constraints.every((g) => g <= 0);

  @override
  String toString() =>
      'Loss(obj=${objective.toStringAsFixed(4)}, '
      'penalty=${penalty.toStringAsFixed(4)}, '
      'V_F=${flutterSpeed.toStringAsFixed(1)}m/s)';
}
