import 'dart:math';

typedef ObjectiveFunction = double Function(List<double> x);

/// Single optimization iteration result.
class OptimizationIteration {
  final int iteration;
  final double bestValue;
  final List<double> bestPoint;
  final bool converged;

  const OptimizationIteration({
    required this.iteration,
    required this.bestValue,
    required this.bestPoint,
    required this.converged,
  });
}

/// Nelder-Mead simplex optimization algorithm.
///
/// Suitable for gradient-free optimization of smooth functions
/// with up to ~20 design variables.
class NelderMead {
  final int maxIterations;
  final double convergenceTolerance;
  final double alpha;   // reflection coefficient (1.0)
  final double gamma;   // expansion coefficient (2.0)
  final double rho;     // contraction coefficient (0.5)
  final double sigma;   // shrink coefficient (0.5)

  const NelderMead({
    this.maxIterations = 500,
    this.convergenceTolerance = 1e-6,
    this.alpha = 1.0,
    this.gamma = 2.0,
    this.rho = 0.5,
    this.sigma = 0.5,
  });

  /// Optimize objective function.
  ///
  /// [objective] - function to minimize
  /// [initialPoint] - starting point (normalized [0,1])
  /// [initialStepSize] - initial simplex size
  ///
  /// Yields [OptimizationIteration] for each iteration (use as stream-like).
  Iterable<OptimizationIteration> optimize({
    required ObjectiveFunction objective,
    required List<double> initialPoint,
    double initialStepSize = 0.1,
  }) sync* {
    final n = initialPoint.length;

    // Build initial simplex (n+1 vertices)
    final simplex = _buildInitialSimplex(initialPoint, initialStepSize);
    var values = simplex.map(objective).toList();

    for (var iter = 0; iter < maxIterations; iter++) {
      // Sort by function value (ascending)
      final order = List.generate(n + 1, (i) => i)
        ..sort((a, b) => values[a].compareTo(values[b]));
      simplex.setAll(0, [for (final i in order) simplex[i]]);
      values = [for (final i in order) values[i]];

      // Check convergence
      final stdDev = _stdDev(values);
      if (stdDev < convergenceTolerance) {
        yield OptimizationIteration(
          iteration: iter,
          bestValue: values[0],
          bestPoint: simplex[0],
          converged: true,
        );
        return;
      }

      yield OptimizationIteration(
        iteration: iter,
        bestValue: values[0],
        bestPoint: simplex[0],
        converged: false,
      );

      // Centroid of all vertices except worst
      final centroid = _centroid(simplex, n);

      // Reflection
      final reflected = _add(centroid, _scale(_sub(centroid, simplex[n]), alpha));
      _clip(reflected);
      final fReflected = objective(reflected);

      if (fReflected < values[0]) {
        // Try expansion
        final expanded = _add(centroid, _scale(_sub(reflected, centroid), gamma));
        _clip(expanded);
        final fExpanded = objective(expanded);
        if (fExpanded < fReflected) {
          simplex[n] = expanded;
          values[n] = fExpanded;
        } else {
          simplex[n] = reflected;
          values[n] = fReflected;
        }
      } else if (fReflected < values[n - 1]) {
        simplex[n] = reflected;
        values[n] = fReflected;
      } else {
        // Contraction
        if (fReflected < values[n]) {
          // Outside contraction
          final contracted = _add(centroid, _scale(_sub(reflected, centroid), rho));
          _clip(contracted);
          final fContracted = objective(contracted);
          if (fContracted <= fReflected) {
            simplex[n] = contracted;
            values[n] = fContracted;
          } else {
            _shrink(simplex, values, objective, n);
          }
        } else {
          // Inside contraction
          final contracted = _add(centroid, _scale(_sub(simplex[n], centroid), rho));
          _clip(contracted);
          final fContracted = objective(contracted);
          if (fContracted < values[n]) {
            simplex[n] = contracted;
            values[n] = fContracted;
          } else {
            _shrink(simplex, values, objective, n);
          }
        }
      }
    }

    // Final result without convergence
    yield OptimizationIteration(
      iteration: maxIterations - 1,
      bestValue: values[0],
      bestPoint: simplex[0],
      converged: false,
    );
  }

  List<List<double>> _buildInitialSimplex(List<double> x0, double step) {
    final n = x0.length;
    final simplex = [List<double>.from(x0)];
    for (var i = 0; i < n; i++) {
      final v = List<double>.from(x0);
      v[i] = (v[i] + step).clamp(0.0, 1.0);
      simplex.add(v);
    }
    return simplex;
  }

  List<double> _centroid(List<List<double>> simplex, int n) {
    final c = List<double>.filled(n, 0.0);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        c[j] += simplex[i][j] / n;
      }
    }
    return c;
  }

  List<double> _add(List<double> a, List<double> b) =>
      List.generate(a.length, (i) => a[i] + b[i]);

  List<double> _sub(List<double> a, List<double> b) =>
      List.generate(a.length, (i) => a[i] - b[i]);

  List<double> _scale(List<double> v, double s) =>
      List.generate(v.length, (i) => v[i] * s);

  void _clip(List<double> v) {
    for (var i = 0; i < v.length; i++) {
      v[i] = v[i].clamp(0.0, 1.0);
    }
  }

  void _shrink(
      List<List<double>> simplex, List<double> values, ObjectiveFunction f, int n) {
    for (var i = 1; i <= n; i++) {
      for (var j = 0; j < simplex[i].length; j++) {
        simplex[i][j] = (simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j]))
            .clamp(0.0, 1.0);
      }
      values[i] = f(simplex[i]);
    }
  }

  double _stdDev(List<double> values) {
    final mean = values.fold(0.0, (s, v) => s + v) / values.length;
    final variance = values.fold(0.0, (s, v) => s + (v - mean) * (v - mean)) / values.length;
    return sqrt(variance);
  }
}
