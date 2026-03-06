import 'dart:math';
import '../../math/matrix.dart';

/// Prandtl-Glauert compressibility correction for subsonic flow.
class PrandtlGlauertCorrection {
  static const double _machLimit = 0.85;

  /// Compressibility factor β = sqrt(1 - M²).
  static double beta(double mach) {
    final m = mach.clamp(0.0, _machLimit);
    return sqrt(1.0 - m * m);
  }

  /// Correct incompressible pressure coefficient: Cp_c = Cp_i / β.
  static double correctCp(double cpIncompressible, double mach) {
    final b = beta(mach);
    if (b < 1e-6) return cpIncompressible;
    return cpIncompressible / b;
  }

  /// Apply P-G correction to an AIC matrix.
  static Matrix correctAIC(Matrix aicIncompressible, double mach) {
    final b = beta(mach);
    if (b < 1e-6) return aicIncompressible;
    return aicIncompressible.scale(1.0 / b);
  }

  /// Check if Mach number is in valid range.
  static bool isValid(double mach) => mach < _machLimit;

  /// Warning message for out-of-range Mach.
  static String? validityWarning(double mach) {
    if (mach >= _machLimit) {
      return 'Prandtl-Glauert correction invalid for M >= ${_machLimit.toStringAsFixed(2)}. '
          'Transonic/supersonic flow requires more advanced corrections.';
    }
    return null;
  }
}
