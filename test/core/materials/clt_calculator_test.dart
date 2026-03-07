import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/materials/clt_calculator.dart';
import 'package:fin_flutter/core/materials/orthotropic_ply.dart';
import 'package:fin_flutter/core/materials/material_database.dart';

void main() {
  group('CLTCalculator', () {
    final calc = CLTCalculator();
    final ply = MaterialDatabase.as43501;

    test('unidirectional laminate: A11 ≈ E1*h, B=0, D11 = E1*h^3/12', () {
      const nPlies = 8;
      final stack = List.generate(nPlies, (_) => LaminatePly(ply, 0.0));
      final abd = calc.compute(stack);

      final h = ply.t * nPlies;
      final E1h = ply.E1 * h; // This is Q11 * h (approximately E1*h for nu12 small)
      final denom = 1.0 - ply.nu12 * ply.nu21;
      final Q11 = ply.E1 / denom;

      expect(abd.A.get(0, 0), closeTo(Q11 * h, Q11 * h * 0.01));
      expect(abd.D.get(0, 0), closeTo(Q11 * h * h * h / 12.0, Q11 * h * h * h / 12.0 * 0.01));

      // B should be zero for symmetric (even number) unidirectional laminate
      // Actually for unidirectional: B = sum of (Q * (z_k^2 - z_{k-1}^2) / 2)
      // For symmetric laminate (even number with midplane symmetry), B = 0.
      // A unidirectional [0]_n is symmetric about midplane:
      for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) {
          // B should be zero for symmetric laminate (even number unidirectional)
          expect(abd.B.get(i, j).abs(), lessThan(1e-3),
              reason: 'B[$i][$j] should be zero for symmetric unidirectional');
        }
      }
    });

    test('symmetric balanced [0/90]_s: A16=A26=0, B=0', () {
      final stack = [
        LaminatePly(ply, 0.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, 0.0),
      ];
      final abd = calc.compute(stack);

      // A16 = A26 = 0 for balanced
      expect(abd.A.get(0, 2).abs(), lessThan(1e-3));
      expect(abd.A.get(1, 2).abs(), lessThan(1e-3));

      // B = 0 for symmetric
      for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) {
          expect(abd.B.get(i, j).abs(), lessThan(1e-3));
        }
      }
    });

    test('cross-ply [0/90]: A11 = A22 (approximately)', () {
      final stack = [
        LaminatePly(ply, 0.0),
        LaminatePly(ply, 90.0),
      ];
      final abd = calc.compute(stack);

      // For a cross-ply with 0 and 90, A11 should equal A22
      expect(abd.A.get(0, 0), closeTo(abd.A.get(1, 1), 1e-3));
    });

    test('quasi-isotropic [0/45/-45/90]_s: A11 = A22, A16 = A26 ≈ 0', () {
      final stack = [
        LaminatePly(ply, 0.0),
        LaminatePly(ply, 45.0),
        LaminatePly(ply, -45.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, 90.0),
        LaminatePly(ply, -45.0),
        LaminatePly(ply, 45.0),
        LaminatePly(ply, 0.0),
      ];
      final abd = calc.compute(stack);

      // A11 ≈ A22 for quasi-isotropic
      expect(abd.A.get(0, 0), closeTo(abd.A.get(1, 1), abd.A.get(0, 0) * 0.01));
      // A16 ≈ 0
      expect(abd.A.get(0, 2).abs(), lessThan(1e-2));
    });

    test('ply angle 45 degrees: nonzero Q16, Q26', () {
      final Q = ply.transformedStiffness(45.0);
      // Off-axis ply should have nonzero Q16 and Q26
      expect(Q.get(0, 2).abs(), greaterThan(1e6));
      expect(Q.get(1, 2).abs(), greaterThan(1e6));
    });

    test('ply angle 0 degrees: Q16=Q26=0', () {
      final Q = ply.transformedStiffness(0.0);
      expect(Q.get(0, 2).abs(), lessThan(1e-3));
      expect(Q.get(1, 2).abs(), lessThan(1e-3));
    });

    test('effective Ex from CLT matches E1 for [0] laminate', () {
      final stack = List.generate(4, (_) => LaminatePly(ply, 0.0));
      final abd = calc.compute(stack);
      final eng = calc.engineeringConstants(abd);
      // For pure [0] laminate, Ex ≈ E1
      expect(eng.Ex, closeTo(ply.E1, ply.E1 * 0.05));
    });
  });
}
