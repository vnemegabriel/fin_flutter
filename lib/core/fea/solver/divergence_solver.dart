import 'dart:math';
import '../../math/matrix.dart';

/// Static aeroelastic divergence solver.
///
/// Divergence occurs when the aerodynamic moment about the root exceeds
/// the structural restoring moment: det([K] - q_D * [A_s]) = 0.
///
/// This is equivalent to finding the minimum positive eigenvalue of
/// [K]^{-1} * [A_s].
class DivergenceSolver {
  /// Solve for divergence dynamic pressure q_D and speed V_D.
  ///
  /// [K] - structural stiffness matrix (n × n)
  /// [A_static] - static aerodynamic stiffness matrix (aerodynamic loads / deflection)
  /// [rho] - air density (kg/m³)
  ///
  /// Returns divergence result with q_D and V_D.
  DivergenceResult solve({
    required Matrix K,
    required Matrix aStatic,
    required double rho,
  }) {
    final n = K.rows;

    // Form B = K^{-1} * A_static
    // Divergence: q_D = minimum positive eigenvalue of B
    Matrix B;
    try {
      B = K.solveMatrix(aStatic);
    } catch (e) {
      return DivergenceResult(
        qDivergence: null,
        vDivergence: null,
        error: 'Singular stiffness matrix: $e',
      );
    }

    // Find minimum positive eigenvalue of B using power iteration on B^{-1}
    // (largest eigenvalue of B^{-1} = 1 / smallest eigenvalue of B)
    final minEigen = _minimumPositiveEigenvalue(B);

    if (minEigen == null || minEigen <= 0) {
      return DivergenceResult(
        qDivergence: null,
        vDivergence: null,
        error: 'No positive divergence eigenvalue found',
      );
    }

    final qD = minEigen;
    final vD = sqrt(2.0 * qD / rho);

    return DivergenceResult(
      qDivergence: qD,
      vDivergence: vD,
    );
  }

  /// Find minimum positive eigenvalue of matrix B using inverse iteration.
  double? _minimumPositiveEigenvalue(Matrix B) {
    final n = B.rows;
    // Use power iteration on B to get the dominant eigenvalue,
    // then use shift-and-invert to get the smallest.

    // Start: find all eigenvalues via QR iteration (small n only)
    if (n > 50) {
      // For large n, use Gershgorin bound + inverse power iteration
      return _inverseIteration(B);
    }

    // For small n, use direct eigenvalue sweep
    return _inverseIteration(B);
  }

  double? _inverseIteration(Matrix B) {
    final n = B.rows;
    const maxIter = 500;
    const tol = 1e-8;

    // Start with shift near 0+ to capture smallest positive eigenvalue
    // Use Rayleigh quotient iteration
    var v = List<double>.generate(n, (i) => i == 0 ? 1.0 : 0.0);
    _normalize(v);

    double lambda = 0.0;
    for (var iter = 0; iter < maxIter; iter++) {
      // v_new = B * v (power iteration for dominant eigenvalue)
      // For minimum eigenvalue, use B^{-1} (shift=0)
      List<double> vNew;
      try {
        vNew = B.solve(v);
      } catch (e) {
        return null;
      }

      final norm = _norm(vNew);
      if (norm < 1e-14) return null;
      _normalize(vNew);

      // Rayleigh quotient: lambda = v^T B v / v^T v
      final Bv = B.multiplyVector(vNew);
      final lambdaNew = _dot(vNew, Bv);

      if ((lambdaNew - lambda).abs() < tol * (1.0 + lambdaNew.abs())) {
        return lambdaNew > 0 ? lambdaNew : null;
      }
      lambda = lambdaNew;
      v = vNew;
    }
    return lambda > 0 ? lambda : null;
  }

  void _normalize(List<double> v) {
    final n = _norm(v);
    if (n > 1e-14) {
      for (var i = 0; i < v.length; i++) v[i] /= n;
    }
  }

  double _norm(List<double> v) {
    double s = 0;
    for (final x in v) s += x * x;
    return sqrt(s);
  }

  double _dot(List<double> a, List<double> b) {
    double s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
}

/// Result of divergence analysis.
class DivergenceResult {
  /// Dynamic pressure at divergence (Pa), or null if not found.
  final double? qDivergence;

  /// Divergence speed V_D (m/s), or null if not found.
  final double? vDivergence;

  /// Error message, if any.
  final String? error;

  const DivergenceResult({
    required this.qDivergence,
    required this.vDivergence,
    this.error,
  });

  bool get hasDivergence => vDivergence != null;

  @override
  String toString() {
    if (hasDivergence) {
      return 'DivergenceResult(V_D=${vDivergence!.toStringAsFixed(1)}m/s, '
          'q_D=${qDivergence!.toStringAsFixed(0)}Pa)';
    }
    return 'DivergenceResult(${error ?? 'no divergence found'})';
  }
}
