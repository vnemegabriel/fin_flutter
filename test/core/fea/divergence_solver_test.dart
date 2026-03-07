import 'dart:math';
import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/fea/solver/divergence_solver.dart';
import 'package:fin_flutter/core/math/matrix.dart';

Matrix _diag(List<double> values) {
  final n = values.length;
  final m = Matrix(n, n);
  for (var i = 0; i < n; i++) m.set(i, i, values[i]);
  return m;
}

void main() {
  group('DivergenceSolver', () {
    final solver = DivergenceSolver();

    test('2×2 diagonal: q_D = minimum of k_i/a_i', () {
      // K = diag(1000, 2000) N/m
      // A_s = diag(2, 1)    N/m per unit q
      // λ₁ = 2/1000 = 0.002  → q_D = 500 Pa
      // λ₂ = 1/2000 = 0.0005 → q_D = 2000 Pa
      // minimum positive: q_D = 500 Pa
      final K = _diag([1000.0, 2000.0]);
      final As = _diag([2.0, 1.0]);
      final rho = 1.225;

      final result = solver.solve(K: K, aStatic: As, rho: rho);

      expect(result.hasDivergence, isTrue);
      expect(result.qDivergence!, closeTo(500.0, 10.0)); // ±2%
      final vExpected = sqrt(2.0 * 500.0 / rho);
      expect(result.vDivergence!, closeTo(vExpected, 1.0));
    });

    test('3×3 system: q_D matches minimum K⁻¹·A_s eigenvalue', () {
      // Use a simple system where K = 3*I, A_s = I
      // K⁻¹ A_s = (1/3) I → all eigenvalues = 1/3
      // q_D = 1/3
      final K = _diag([3.0, 3.0, 3.0]);
      final As = _diag([1.0, 1.0, 1.0]);
      final rho = 1.225;

      final result = solver.solve(K: K, aStatic: As, rho: rho);

      expect(result.hasDivergence, isTrue);
      expect(result.qDivergence!, closeTo(1.0 / 3.0, 0.02));
    });

    test('q_D is strictly positive for physical inputs', () {
      // Any physical system with K > 0 and A_s > 0 must yield q_D > 0
      final K = _diag([500.0, 800.0, 1200.0]);
      final As = _diag([1.0, 0.5, 0.25]);
      final rho = 1.0;

      final result = solver.solve(K: K, aStatic: As, rho: rho);

      if (result.hasDivergence) {
        expect(result.qDivergence!, greaterThan(0.0));
        expect(result.vDivergence!, greaterThan(0.0));
      }
    });

    test('returns no-divergence gracefully when A_s is negative-definite', () {
      // A_s = diag(−1, −1) → K⁻¹ A_s has negative eigenvalues → no divergence
      final K = _diag([1000.0, 1000.0]);
      final As = _diag([-1.0, -1.0]);
      final rho = 1.225;

      final result = solver.solve(K: K, aStatic: As, rho: rho);

      // Either no divergence found or hasDivergence is false — no crash
      expect(result.hasDivergence, isFalse);
      expect(result.qDivergence, isNull);
      expect(result.error, isNotNull);
    });

    test('DivergenceResult toString is informative', () {
      final resultWithDiv = DivergenceResult(
        qDivergence: 500.0,
        vDivergence: 28.57,
      );
      expect(resultWithDiv.toString(), contains('28.6'));

      final resultNoDiv = DivergenceResult(
        qDivergence: null,
        vDivergence: null,
        error: 'No positive eigenvalue',
      );
      expect(resultNoDiv.toString(), contains('No positive eigenvalue'));
    });
  });
}
