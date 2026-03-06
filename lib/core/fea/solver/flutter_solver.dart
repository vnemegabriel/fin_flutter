import 'dart:math';
import '../../math/matrix.dart';
import '../../math/complex_number.dart';
import '../../math/complex_matrix.dart';

/// A single point on the V-g (or V-f) diagram.
class VGPoint {
  final double velocity; // m/s
  final double frequency; // Hz
  final double damping;  // structural damping g (dimensionless)
  final int modeIndex;

  const VGPoint({
    required this.velocity,
    required this.frequency,
    required this.damping,
    required this.modeIndex,
  });
}

/// Flutter analysis result.
class FlutterResult {
  /// Flutter speed (m/s), or null if no flutter found in range.
  final double? flutterSpeed;

  /// Divergence speed (m/s) from static aeroelastic analysis.
  final double? divergenceSpeed;

  /// Flutter mode index.
  final int? flutterMode;

  /// Flutter frequency (Hz) at flutter point.
  final double? flutterFrequency;

  /// V-g curves for all modes.
  final List<List<VGPoint>> vgCurves;

  /// V-f curves for all modes.
  final List<List<VGPoint>> vfCurves;

  const FlutterResult({
    this.flutterSpeed,
    this.divergenceSpeed,
    this.flutterMode,
    this.flutterFrequency,
    required this.vgCurves,
    required this.vfCurves,
  });

  bool get hasFlutter => flutterSpeed != null;
  bool get hasDivergence => divergenceSpeed != null;

  @override
  String toString() {
    if (hasFlutter) {
      return 'FlutterResult(V_F=${flutterSpeed!.toStringAsFixed(1)}m/s, '
          'f=${flutterFrequency!.toStringAsFixed(2)}Hz)';
    }
    return 'FlutterResult(no flutter in range)';
  }
}

/// Aeroelastic flutter solver implementing the U-g method.
///
/// The U-g method sweeps dynamic pressure (velocity) and finds where
/// structural damping g crosses zero, indicating the flutter boundary.
class FlutterSolver {
  /// Solve for flutter using the U-g (k-method) approach.
  ///
  /// [K_modal] - modal stiffness (diagonal, n_modes × n_modes)
  /// [M_modal] - modal mass (diagonal, n_modes × n_modes)
  /// [Q_modal] - generalized aerodynamic force matrix (n_modes × n_modes)
  ///             Q_modal = Phi^T * [aerodynamic_stiffness] * Phi
  /// [rho] - air density (kg/m³)
  /// [velocities] - velocity sweep values (m/s)
  FlutterResult solveUG({
    required Matrix kModal,
    required Matrix mModal,
    required Matrix qModal,
    required double rho,
    required List<double> velocities,
  }) {
    final nModes = kModal.rows;
    final nVels = velocities.length;

    // VG curves: for each mode, list of (V, g, freq) points
    final vgCurves = List.generate(nModes, (_) => <VGPoint>[]);
    final vfCurves = List.generate(nModes, (_) => <VGPoint>[]);

    // Track g values for flutter crossing detection
    final prevDamping = List<double?>.filled(nModes, null);
    double? flutterSpeed;
    double? flutterFrequency;
    int? flutterMode;

    for (final V in velocities) {
      final q = 0.5 * rho * V * V; // dynamic pressure

      // Flutter equation (U-g formulation):
      // [K_g + i*g*K] - omega^2*[M] - q*[Q] = 0
      // Rearranged as eigenvalue problem in g:
      // Let z = 1/omega^2, then:
      // ([K] - q*[Q]) * z - [M] = -i*g*[K]*z (artificial damping)
      //
      // For the k-method, we assume harmonic motion at some k=omega*b/V
      // and find g needed for equilibrium.
      //
      // Simplified approach: form the aeroelastic eigenvalue matrix
      // and extract complex eigenvalues.

      // Form: A = [K] - q*[Q], solve A*{q} = omega^2*[M]*{q}
      // i.e., [M]^{-1}*([K] - q*[Q]) * {phi} = omega^2 * {phi}
      // But this ignores the imaginary part of Q (unsteady aerodynamics).
      // For steady aerodynamic forces (quasi-steady), Q_modal is real.

      // Form system matrix: Z = M^{-1} * (K - q*Q)
      final Z = Matrix(nModes, nModes);
      for (var i = 0; i < nModes; i++) {
        final mi = mModal.get(i, i);
        if (mi.abs() < 1e-14) continue;
        for (var j = 0; j < nModes; j++) {
          Z.set(i, j, (kModal.get(i, j) - q * qModal.get(i, j)) / mi);
        }
      }

      // Extract eigenvalues of Z (complex in general due to off-diagonal terms)
      // For simplified quasi-steady: eigenvalues are omega^2 (approximately real)
      final evals = _realEigenvalues(Z, nModes);

      for (var m = 0; m < nModes; m++) {
        final lambda = evals[m];
        double freq, damping;

        if (lambda >= 0) {
          // Stable oscillatory
          freq = sqrt(lambda) / (2.0 * pi);
          damping = 0.0; // quasi-steady: no aerodynamic damping (conservative)
        } else {
          // Aeroelastic divergence-like: real eigenvalue negative
          freq = 0.0;
          damping = 1.0; // artificially large damping for display
        }

        // For U-g: estimate damping from imaginary part of eigenvalue
        // This requires complex arithmetic for off-diagonal Q. We estimate:
        // g ≈ -Im(omega^2) / |omega^2| (for lightly damped modes)
        // With quasi-steady assumption, use flutter margin indicator:
        damping = _estimateDamping(kModal, mModal, qModal, m, q);

        vgCurves[m].add(VGPoint(
          velocity: V,
          frequency: freq,
          damping: damping,
          modeIndex: m,
        ));
        vfCurves[m].add(VGPoint(
          velocity: V,
          frequency: freq,
          damping: damping,
          modeIndex: m,
        ));

        // Detect flutter crossing (g goes from negative to positive)
        final prevG = prevDamping[m];
        if (prevG != null && prevG < 0 && damping >= 0 && flutterSpeed == null) {
          // Linear interpolation to find crossing velocity
          final vPrev = velocities[velocities.indexOf(V) - 1];
          final vFlutter = vPrev + (V - vPrev) * (-prevG) / (damping - prevG);
          flutterSpeed = vFlutter;
          flutterMode = m;
          flutterFrequency = freq;
        }
        prevDamping[m] = damping;
      }
    }

    return FlutterResult(
      flutterSpeed: flutterSpeed,
      flutterMode: flutterMode,
      flutterFrequency: flutterFrequency,
      vgCurves: vgCurves,
      vfCurves: vfCurves,
    );
  }

  /// Estimate structural damping using the generalized aerodynamic force approach.
  ///
  /// For a single-mode quasi-steady analysis:
  /// omega^2 = (k_ii - q * Q_ii) / m_ii
  /// The "aerodynamic damping" indicator:
  /// g ≈ -(q * Q_ii) / k_ii + 1 - 1 = q*Q_ii/k_ii - 1 (flutter when g >= 0)
  ///
  /// This is the ratio of aerodynamic stiffness to structural stiffness.
  double _estimateDamping(
      Matrix K, Matrix M, Matrix Q, int modeIdx, double q) {
    final kii = K.get(modeIdx, modeIdx);
    final qii = Q.get(modeIdx, modeIdx);
    if (kii.abs() < 1e-14) return 0.0;
    // g < 0: stable, g > 0: flutter
    return q * qii / kii - 1.0;
  }

  /// Compute approximate real eigenvalues of matrix Z.
  List<double> _realEigenvalues(Matrix Z, int n) {
    // For a symmetric matrix, use power iteration / Gershgorin-based approach.
    // For small n, use direct computation via characteristic polynomial.
    // Here we use the diagonal approximation for quasi-steady flutter:
    return List.generate(n, (i) => Z.get(i, i));
  }
}
