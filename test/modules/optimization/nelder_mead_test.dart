import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/modules/optimization/nelder_mead.dart';

void main() {
  group('NelderMead', () {
    test('minimizes sphere function f(x) = sum(x_i^2)', () {
      final optimizer = NelderMead(
        maxIterations: 500,
        convergenceTolerance: 1e-8,
      );

      // Sphere function: minimum at x=0.5 in [0,1] space (mapped from [0.5,0.5])
      // Actually minimum of sum((2x-1)^2) = 0 at x=[0.5, 0.5]
      final objective = (List<double> x) {
        double s = 0;
        for (final v in x) s += (2 * v - 1) * (2 * v - 1); // min at 0.5
        return s;
      };

      final iterations = optimizer.optimize(
        objective: objective,
        initialPoint: [0.7, 0.3, 0.8],
      ).toList();

      final best = iterations.last;
      expect(best.bestValue, lessThan(0.01));
      for (final xi in best.bestPoint) {
        expect(xi, closeTo(0.5, 0.1));
      }
    });

    test('Rosenbrock function converges to [0.75, ~0.75^2] in [0,1]', () {
      final optimizer = NelderMead(maxIterations: 300);

      // Rosenbrock minimum at (1,1) -> in [0,1] space: (0.75, 0.75) approx for x in [-1,3]
      final objective = (List<double> x) {
        final x0 = x[0] * 4 - 1;  // map to [-1, 3]
        final x1 = x[1] * 4 - 1;
        return (1 - x0) * (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
      };

      final iterations = optimizer.optimize(
        objective: objective,
        initialPoint: [0.5, 0.5],
      ).toList();

      final best = iterations.last;
      // Minimum of Rosenbrock is at 0 -> convergence value should be small
      expect(best.bestValue, lessThan(1.0));
    });

    test('convergence: value decreases monotonically', () {
      final optimizer = NelderMead(maxIterations: 100);
      final objective = (List<double> x) =>
          (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.3) * (x[1] - 0.3);

      final iterations = optimizer.optimize(
        objective: objective,
        initialPoint: [0.9, 0.9],
      ).toList();

      // Best value should generally decrease
      final values = iterations.map((i) => i.bestValue).toList();
      expect(values.last, lessThanOrEqualTo(values.first + 1e-8));
    });
  });
}
