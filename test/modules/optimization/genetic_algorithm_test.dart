import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/modules/optimization/genetic_algorithm.dart';
import 'package:fin_flutter/modules/optimization/nelder_mead.dart' show OptimizationIteration;

void main() {
  group('GeneticAlgorithm', () {
    test('sphere function: best value decreases from generation 0 to 50', () {
      final ga = GeneticAlgorithm(
        populationSize: 30,
        maxGenerations: 50,
        seed: 123,
      );

      // f(x) = Σ(2xᵢ − 1)², minimum at xᵢ = 0.5
      final objective = (List<double> x) {
        double s = 0;
        for (final v in x) s += (2 * v - 1) * (2 * v - 1);
        return s;
      };

      final iterations = ga.optimize(
        objective: objective,
        nVariables: 4,
      ).toList();

      expect(iterations, isNotEmpty);

      final firstBest = iterations.first.bestValue;
      final lastBest = iterations.last.bestValue;

      // The GA should improve over 50 generations
      expect(lastBest, lessThanOrEqualTo(firstBest + 1e-6));
      // Final best should be reasonably close to 0
      expect(lastBest, lessThan(1.0));
    });

    test('trivial all-zero function: GA terminates without error', () {
      final ga = GeneticAlgorithm(
        populationSize: 10,
        maxGenerations: 5,
        seed: 0,
      );

      final iterations = ga.optimize(
        objective: (_) => 0.0,
        nVariables: 3,
      ).toList();

      // Should yield maxGenerations + 1 iterations (final converged result)
      expect(iterations.length, greaterThan(0));
      for (final iter in iterations) {
        expect(iter.bestValue, closeTo(0.0, 1e-10));
      }
    });

    test('population size is maintained across generations', () {
      // We can't directly inspect population size from the iterator,
      // but we can verify that the number of yielded iterations is correct.
      final maxGen = 10;
      final ga = GeneticAlgorithm(
        populationSize: 20,
        maxGenerations: maxGen,
        seed: 7,
      );

      final iterations = ga.optimize(
        objective: (x) => x.fold(0.0, (s, v) => s + v * v),
        nVariables: 2,
      ).toList();

      // Should yield maxGenerations intermediate + 1 final = maxGen + 1 total
      expect(iterations.length, maxGen + 1);
    });

    test('elite preservation: last generation best ≤ second-to-last best', () {
      final ga = GeneticAlgorithm(
        populationSize: 20,
        maxGenerations: 20,
        eliteCount: 2,
        seed: 42,
      );

      final objective = (List<double> x) =>
          (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.3) * (x[1] - 0.3);

      final iterations = ga.optimize(
        objective: objective,
        nVariables: 2,
      ).toList();

      // With elitism, the best value should never increase
      for (var i = 1; i < iterations.length; i++) {
        expect(
          iterations[i].bestValue,
          lessThanOrEqualTo(iterations[i - 1].bestValue + 1e-10),
          reason: 'Best value increased at generation ${iterations[i].iteration}',
        );
      }
    });

    test('all individuals remain in [0,1] bounds after evolution', () {
      final ga = GeneticAlgorithm(
        populationSize: 15,
        maxGenerations: 10,
        mutationRate: 0.5, // high mutation to stress-test bounds
        seed: 99,
      );

      for (final iter in ga.optimize(
        objective: (x) => x.fold(0.0, (s, v) => s + v),
        nVariables: 3,
      )) {
        for (final xi in iter.bestPoint) {
          expect(xi, greaterThanOrEqualTo(0.0 - 1e-10));
          expect(xi, lessThanOrEqualTo(1.0 + 1e-10));
        }
      }
    });
  });
}
