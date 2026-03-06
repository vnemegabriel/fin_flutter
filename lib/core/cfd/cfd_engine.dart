import '../math/matrix.dart';
import '../models/fin_geometry.dart';
import '../models/flight_condition.dart';
import 'corrections/prandtl_glauert.dart';
import 'vlm/vortex_lattice.dart';

/// Configuration for CFD/VLM analysis.
class CFDConfig {
  final int chordwisePanels;
  final int spanwisePanels;
  final bool applyCompressibilityCorrection;
  final double angleOfAttackDeg;

  const CFDConfig({
    this.chordwisePanels = 8,
    this.spanwisePanels = 12,
    this.applyCompressibilityCorrection = true,
    this.angleOfAttackDeg = 2.0,
  });
}

/// Orchestrates CFD analysis pipeline.
class CFDEngine {
  final CFDConfig config;

  const CFDEngine({this.config = const CFDConfig()});

  /// Run VLM analysis and return aerodynamic results + AIC matrix.
  CFDResult analyze({
    required FinGeometry geometry,
    required FlightCondition condition,
  }) {
    final vlm = VortexLattice();
    final aoa = config.angleOfAttackDeg * 3.14159265358979 / 180.0;

    VLMResult vlmResult = vlm.solve(
      geometry: geometry,
      condition: condition,
      chordwisePanels: config.chordwisePanels,
      spanwisePanels: config.spanwisePanels,
      angleOfAttackRad: aoa,
    );

    // Apply compressibility correction
    Matrix aicCorrected = vlmResult.aicMatrix;
    if (config.applyCompressibilityCorrection &&
        PrandtlGlauertCorrection.isValid(condition.mach)) {
      aicCorrected =
          PrandtlGlauertCorrection.correctAIC(vlmResult.aicMatrix, condition.mach);
    }

    return CFDResult(
      vlmResult: vlmResult,
      aicCorrected: aicCorrected,
      compressibilityCorrectionApplied:
          config.applyCompressibilityCorrection &&
              PrandtlGlauertCorrection.isValid(condition.mach),
      validityWarning: PrandtlGlauertCorrection.validityWarning(condition.mach),
    );
  }
}

/// Results from CFD analysis.
class CFDResult {
  final VLMResult vlmResult;
  final Matrix aicCorrected;
  final bool compressibilityCorrectionApplied;
  final String? validityWarning;

  const CFDResult({
    required this.vlmResult,
    required this.aicCorrected,
    required this.compressibilityCorrectionApplied,
    this.validityWarning,
  });

  double get CL => vlmResult.CL;
  double get CDi => vlmResult.CDi;
}
