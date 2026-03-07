import 'design_variables.dart';
import 'nelder_mead.dart';

/// Result of an optimization run.
class OptimizationResult {
  final DesignVariables optimalVariables;
  final double optimalObjective;
  final double optimalFlutterSpeed;
  final double optimalDivergenceSpeed;
  final double optimalMass;
  final List<OptimizationIteration> history;
  final bool converged;
  final int totalIterations;
  final String optimizerName;

  const OptimizationResult({
    required this.optimalVariables,
    required this.optimalObjective,
    required this.optimalFlutterSpeed,
    required this.optimalDivergenceSpeed,
    required this.optimalMass,
    required this.history,
    required this.converged,
    required this.totalIterations,
    required this.optimizerName,
  });

  @override
  String toString() =>
      'OptimizationResult(obj=${optimalObjective.toStringAsFixed(4)}, '
      'V_F=${optimalFlutterSpeed.toStringAsFixed(1)}m/s, '
      'converged=$converged after $totalIterations iters)';
}
