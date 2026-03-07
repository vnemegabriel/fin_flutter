import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/models/flight_condition.dart';
import 'package:fin_flutter/modules/optimization/loss_function.dart';

void main() {
  group('LossEvaluation', () {
    test('feasible is true when all constraints are satisfied (≤ 0)', () {
      final eval = LossEvaluation(
        objective: 0.5,
        penalty: 0.0,
        total: 0.5,
        flutterSpeed: 250.0,
        divergenceSpeed: 400.0,
        mass: 0.05,
        constraints: [-0.1, -0.5, -1.0],
      );

      expect(eval.feasible, isTrue);
    });

    test('feasible is false when any constraint is violated (> 0)', () {
      final eval = LossEvaluation(
        objective: 0.5,
        penalty: 10.0,
        total: 10.5,
        flutterSpeed: 80.0,
        divergenceSpeed: 200.0,
        mass: 0.05,
        constraints: [-0.1, 0.5, -1.0], // constraint[1] violated
      );

      expect(eval.feasible, isFalse);
    });

    test('total equals objective + penalty', () {
      final eval = LossEvaluation(
        objective: 1.25,
        penalty: 3.75,
        total: 5.0,
        flutterSpeed: 100.0,
        divergenceSpeed: 200.0,
        mass: 0.1,
        constraints: [0.0],
      );

      expect(eval.total, closeTo(eval.objective + eval.penalty, 1e-12));
    });

    test('flutter loss term is zero when V_F equals V_ref', () {
      // w1*(V_F/V_ref − 1)² = 0 when V_F = V_ref
      const vRef = 200.0;
      const vFlutter = vRef;
      final flutterTerm = (vFlutter / vRef - 1.0) * (vFlutter / vRef - 1.0);

      expect(flutterTerm, closeTo(0.0, 1e-15));
    });

    test('penalty is zero and feasible when all constraints are zero', () {
      final eval = LossEvaluation(
        objective: 0.0,
        penalty: 0.0,
        total: 0.0,
        flutterSpeed: 200.0,
        divergenceSpeed: 200.0,
        mass: 1.0,
        constraints: [0.0, 0.0, 0.0],
      );

      expect(eval.penalty, closeTo(0.0, 1e-15));
      // g_i = 0 counts as satisfied
      expect(eval.feasible, isTrue);
    });

    test('toString contains flutter speed value', () {
      final eval = LossEvaluation(
        objective: 0.1,
        penalty: 0.0,
        total: 0.1,
        flutterSpeed: 175.3,
        divergenceSpeed: 300.0,
        mass: 0.05,
        constraints: [-1.0],
      );

      expect(eval.toString(), contains('175.3'));
    });
  });

  group('FlutterLossFunction default configuration', () {
    final condition = FlightCondition.isa(altitude: 0, mach: 0.3);

    test('default weights are physically sensible', () {
      final lf = FlutterLossFunction(
        vRef: 200.0,
        flightCondition: condition,
      );

      expect(lf.w1, greaterThan(0.0));
      expect(lf.w2, greaterThan(0.0));
      expect(lf.w3, greaterThan(0.0));
      expect(lf.mu, greaterThan(0.0));
      // Flutter objective should be weighted at least as much as mass penalty
      expect(lf.w1, greaterThanOrEqualTo(lf.w3));
    });

    test('penalty coefficient mu is positive', () {
      final lf = FlutterLossFunction(vRef: 300.0, flightCondition: condition);
      expect(lf.mu, greaterThan(0.0));
    });
  });
}
