import 'dart:math';
import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/fea/solver/eigenvalue_solver.dart';
import 'package:fin_flutter/core/math/matrix.dart';

void main() {
  group('EigenvalueSolver', () {
    final solver = EigenvalueSolver();

    test('2x2 symmetric: eigenvalues correct', () {
      // K = [3 1; 1 3], M = I
      // Eigenvalues: 2, 4
      final K = Matrix.fromList([[3.0, 1.0], [1.0, 3.0]]);
      final M = Matrix.identity(2);

      final result = solver.solveGeneralized(K: K, M: M, nModes: 2);
      final evals = result.eigenvalues..sort();

      expect(evals[0], closeTo(2.0, 0.1));
      expect(evals[1], closeTo(4.0, 0.1));
    });

    test('diagonal K and M: eigenvalues = K_ii / M_ii', () {
      final K = Matrix.diagonal([100.0, 400.0, 900.0]);
      final M = Matrix.diagonal([1.0, 1.0, 1.0]);

      final result = solver.solveGeneralized(K: K, M: M, nModes: 3);
      final evals = result.eigenvalues..sort();

      expect(evals[0], closeTo(100.0, 5.0));
      expect(evals[1], closeTo(400.0, 20.0));
      expect(evals[2], closeTo(900.0, 45.0));
    });

    test('natural frequencies from eigenvalues', () {
      // omega^2 = 1000 -> omega = 31.62 rad/s -> f = 5.03 Hz
      final K = Matrix.diagonal([1000.0, 4000.0]);
      final M = Matrix.diagonal([1.0, 1.0]);

      final result = solver.solveGeneralized(K: K, M: M, nModes: 2);
      final freqs = result.naturalFrequenciesHz;

      expect(freqs.first, closeTo(sqrt(1000) / (2 * pi), 1.0));
    });
  });
}
