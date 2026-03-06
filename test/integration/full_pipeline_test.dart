import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/models/fin_geometry.dart';
import 'package:fin_flutter/core/models/flight_condition.dart';
import 'package:fin_flutter/core/materials/material_database.dart';
import 'package:fin_flutter/core/materials/laminate.dart';
import 'package:fin_flutter/core/fea/fea_engine.dart';
import 'package:fin_flutter/core/fea/mesh/mesh.dart';
import 'package:fin_flutter/core/cfd/cfd_engine.dart';
import 'package:fin_flutter/modules/flutter_analysis/flutter_analysis_module.dart';
import 'package:fin_flutter/modules/flutter_analysis/flutter_input.dart';

void main() {
  group('Full Pipeline Integration', () {
    late FinGeometry geometry;
    late FlightCondition condition;

    setUp(() {
      geometry = FinGeometry(
        span: 0.15,
        rootChord: 0.12,
        tipChord: 0.06,
        sweepLength: 0.04,
        thickness: 0.003,
      );
      condition = FlightCondition.isa(altitude: 1000, mach: 0.3);
    });

    test('FEA produces positive natural frequencies', () {
      final ply = MaterialDatabase.as43501;
      final lam = Laminate.quasiIsotropic(ply);
      final abd = lam.computeABD();

      final engine = FEAEngine(
        config: FEAConfig(
          spanwiseElements: 3,
          chordwiseElements: 2,
          elementType: ElementType.kirchhoffDKQ,
          nModes: 3,
        ),
      );

      final result = engine.analyze(geometry: geometry, abd: abd);

      expect(result.nModes, greaterThan(0));
      // All valid modes should have positive frequencies
      for (final freq in result.naturalFrequenciesHz) {
        expect(freq, greaterThanOrEqualTo(0));
      }
    });

    test('VLM produces finite CL for positive AoA', () {
      final engine = CFDEngine(
        config: const CFDConfig(
          chordwisePanels: 3,
          spanwisePanels: 3,
          angleOfAttackDeg: 3.0,
        ),
      );
      final result = engine.analyze(geometry: geometry, condition: condition);

      expect(result.CL.isFinite, true);
      expect(result.CL, greaterThan(0));
    });

    test('Full flutter analysis completes without error', () {
      final ply = MaterialDatabase.t3005208;
      final lam = Laminate.crossPly(ply, 2);
      final abd = lam.computeABD();

      final module = FlutterAnalysisModule();
      final output = module.analyze(FlutterAnalysisInput(
        geometry: geometry,
        laminate: abd,
        flightCondition: condition,
        feaConfig: const FEAConfig(
          spanwiseElements: 3,
          chordwiseElements: 2,
          nModes: 2,
        ),
        cfdConfig: const CFDConfig(
          chordwisePanels: 2,
          spanwisePanels: 3,
        ),
        velocityMin: 20.0,
        velocityMax: 300.0,
        velocitySteps: 10,
        maxFlightVelocity: 100.0,
      ));

      // Should complete without throwing
      expect(output.feaResult.nModes, greaterThan(0));
      expect(output.cfdResult.CL.isFinite, true);
      // V-g curves should have data
      expect(output.flutterResult.vgCurves.isNotEmpty, true);
    });

    test('Geometry derived properties are physically correct', () {
      expect(geometry.planformArea, greaterThan(0));
      expect(geometry.aspectRatio, greaterThan(0));
      expect(geometry.sweepAngleDeg, greaterThan(0));
      expect(geometry.taperRatio, lessThan(1.0));
      expect(geometry.taperRatio, greaterThan(0.0));
    });

    test('ISA atmosphere gives correct sea level density', () {
      final seaLevel = FlightCondition.isa(altitude: 0, mach: 0.3);
      expect(seaLevel.density, closeTo(1.225, 0.01));
      expect(seaLevel.speedOfSound, closeTo(340.3, 2.0));
    });
  });
}
