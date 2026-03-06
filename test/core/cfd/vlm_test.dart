import 'dart:math';
import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/cfd/vlm/vortex_lattice.dart';
import 'package:fin_flutter/core/models/fin_geometry.dart';
import 'package:fin_flutter/core/models/flight_condition.dart';

void main() {
  group('VortexLattice', () {
    test('lift coefficient positive for positive angle of attack', () {
      final geo = FinGeometry(
        span: 0.2,
        rootChord: 0.15,
        tipChord: 0.08,
        sweepLength: 0.05,
        thickness: 0.003,
      );
      final condition = FlightCondition.isa(altitude: 0, mach: 0.3);
      final vlm = VortexLattice();

      final result = vlm.solve(
        geometry: geo,
        condition: condition,
        chordwisePanels: 4,
        spanwisePanels: 4,
        angleOfAttackRad: 3.0 * pi / 180.0,
      );

      // CL should be positive for positive AoA on a flat plate fin
      expect(result.CL, greaterThan(0));
    });

    test('AIC matrix is square with correct dimension', () {
      final geo = FinGeometry(
        span: 0.15,
        rootChord: 0.12,
        tipChord: 0.06,
        sweepLength: 0.04,
        thickness: 0.003,
      );
      final condition = FlightCondition.isa(altitude: 0, mach: 0.2);
      final vlm = VortexLattice();

      final result = vlm.solve(
        geometry: geo,
        condition: condition,
        chordwisePanels: 3,
        spanwisePanels: 3,
        angleOfAttackRad: 0.0,
      );

      final nPanels = 3 * 3;
      expect(result.aicMatrix.rows, nPanels);
      expect(result.aicMatrix.cols, nPanels);
    });

    test('CL increases with angle of attack', () {
      final geo = FinGeometry(
        span: 0.2,
        rootChord: 0.15,
        tipChord: 0.1,
        sweepLength: 0.0,
        thickness: 0.003,
      );
      final condition = FlightCondition.isa(altitude: 0, mach: 0.2);
      final vlm = VortexLattice();

      final r1 = vlm.solve(
        geometry: geo, condition: condition,
        chordwisePanels: 3, spanwisePanels: 3,
        angleOfAttackRad: 2.0 * pi / 180.0,
      );
      final r2 = vlm.solve(
        geometry: geo, condition: condition,
        chordwisePanels: 3, spanwisePanels: 3,
        angleOfAttackRad: 5.0 * pi / 180.0,
      );

      expect(r2.CL, greaterThan(r1.CL));
    });
  });
}
