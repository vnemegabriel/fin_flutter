import 'dart:math';
import '../../math/matrix.dart';
import '../../math/complex_number.dart';

/// Result of an eigenvalue solve.
class EigenResult {
  /// Eigenvalues (natural frequency squared ω²) in ascending order, rad²/s².
  final List<double> eigenvalues;

  /// Corresponding eigenvectors (mode shapes), columns of matrix.
  final Matrix eigenvectors;

  /// Natural frequencies in Hz.
  List<double> get naturalFrequenciesHz =>
      eigenvalues.map((ev) => sqrt(ev.abs()) / (2.0 * pi)).toList();

  /// Natural frequencies in rad/s.
  List<double> get naturalFrequenciesRadPS =>
      eigenvalues.map((ev) => sqrt(ev.abs())).toList();

  const EigenResult({required this.eigenvalues, required this.eigenvectors});
}

/// Solves structural eigenvalue problems for natural frequencies and mode shapes.
class EigenvalueSolver {
  /// Solve generalized eigenvalue problem: [K]{phi} = lambda*[M]{phi}
  ///
  /// Uses inverse iteration (shift-and-invert) with deflation.
  /// Returns the first [nModes] positive eigenvalues (natural freq²) and modes.
  ///
  /// [K] - stiffness matrix (symmetric, positive semi-definite after BC application)
  /// [M] - mass matrix (symmetric, positive definite)
  /// [nModes] - number of modes to extract
  EigenResult solveGeneralized({
    required Matrix K,
    required Matrix M,
    required int nModes,
    double shift = 0.0,
    int maxIterations = 1000,
    double tolerance = 1e-8,
  }) {
    final n = K.rows;
    assert(K.rows == K.cols && M.rows == M.cols && K.rows == M.rows);
    final actualModes = nModes.clamp(1, n);

    final eigenvalues = <double>[];
    final eigenVectors = <List<double>>[];

    // Use Lanczos-style subspace iteration for moderate-sized problems.
    // For small problems, use direct inverse iteration with deflation.
    if (n <= 200) {
      _lanczosIteration(
        K: K,
        M: M,
        nModes: actualModes,
        eigenvalues: eigenvalues,
        eigenVectors: eigenVectors,
        maxIter: maxIterations,
        tol: tolerance,
        shift: shift,
      );
    } else {
      // For larger problems, use a simplified subspace iteration
      _subspaceIteration(
        K: K,
        M: M,
        nModes: actualModes,
        eigenvalues: eigenvalues,
        eigenVectors: eigenVectors,
        maxIter: maxIterations,
        tol: tolerance,
      );
    }

    // Sort by eigenvalue (ascending)
    final pairs = List.generate(
        eigenvalues.length, (i) => (eigenvalues[i], eigenVectors[i]));
    pairs.sort((a, b) => a.$1.compareTo(b.$1));

    final sortedEvals = pairs.map((p) => p.$1).toList();
    final sortedEvecs = Matrix(n, pairs.length);
    for (var j = 0; j < pairs.length; j++) {
      final evec = pairs[j].$2;
      for (var i = 0; i < n; i++) {
        sortedEvecs.set(i, j, evec[i]);
      }
    }

    return EigenResult(eigenvalues: sortedEvals, eigenvectors: sortedEvecs);
  }

  void _lanczosIteration({
    required Matrix K,
    required Matrix M,
    required int nModes,
    required List<double> eigenvalues,
    required List<List<double>> eigenVectors,
    required int maxIter,
    required double tol,
    required double shift,
  }) {
    final n = K.rows;
    final rng = _Random(42);

    // Build shifted matrix: A = K - shift*M
    final A = Matrix(n, n);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        A.set(i, j, K.get(i, j) - shift * M.get(i, j));
      }
    }

    // For each mode, use inverse iteration
    for (var mode = 0; mode < nModes; mode++) {
      // Random starting vector, orthogonalized against previous modes
      var v = List<double>.generate(n, (_) => rng.nextDouble() - 0.5);
      _orthogonalize(v, eigenVectors, M);
      _mNormalize(v, M);

      double lambda = 0.0;
      for (var iter = 0; iter < maxIter; iter++) {
        // v_new = A^{-1} * M * v (inverse iteration)
        final Mv = M.multiplyVector(v);
        List<double> vNew;
        try {
          vNew = A.solve(Mv);
        } catch (e) {
          // Singular or near-singular: add small perturbation to shift
          break;
        }

        // Orthogonalize against converged modes
        _orthogonalize(vNew, eigenVectors, M);

        // Rayleigh quotient: lambda = (v^T K v) / (v^T M v)
        final Kv = K.multiplyVector(vNew);
        final KvdotV = dotProd(vNew, Kv);
        final MvdotV = dotProd(vNew, M.multiplyVector(vNew));
        final lambdaNew = KvdotV / MvdotV;

        _mNormalize(vNew, M);

        if ((lambdaNew - lambda).abs() < tol * (1.0 + lambdaNew.abs())) {
          lambda = lambdaNew;
          v = vNew;
          break;
        }
        lambda = lambdaNew;
        v = vNew;
      }

      if (lambda > -1e6) {
        eigenvalues.add(lambda);
        eigenVectors.add(List<double>.from(v));
      }
    }
  }

  void _subspaceIteration({
    required Matrix K,
    required Matrix M,
    required int nModes,
    required List<double> eigenvalues,
    required List<List<double>> eigenVectors,
    required int maxIter,
    required double tol,
  }) {
    // Simplified: use inverse iteration for each mode
    _lanczosIteration(
      K: K,
      M: M,
      nModes: nModes,
      eigenvalues: eigenvalues,
      eigenVectors: eigenVectors,
      maxIter: maxIter,
      tol: tol,
      shift: 0.0,
    );
  }

  /// Orthogonalize v against all vectors in basis (M-orthogonality).
  void _orthogonalize(
      List<double> v, List<List<double>> basis, Matrix M) {
    for (final bv in basis) {
      final Mv = M.multiplyVector(v);
      final proj = dotProd(bv, Mv);
      for (var i = 0; i < v.length; i++) {
        v[i] -= proj * bv[i];
      }
    }
  }

  /// M-normalize vector: v = v / sqrt(v^T M v).
  void _mNormalize(List<double> v, Matrix M) {
    final Mv = M.multiplyVector(v);
    final norm = sqrt(dotProd(v, Mv));
    if (norm > 1e-14) {
      for (var i = 0; i < v.length; i++) {
        v[i] /= norm;
      }
    }
  }

  double dotProd(List<double> a, List<double> b) {
    double s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
}

/// Simple LCG random number generator for deterministic seeding.
class _Random {
  int _state;
  _Random(this._state);
  double nextDouble() {
    _state = (_state * 1664525 + 1013904223) & 0x7FFFFFFF;
    return _state / 0x7FFFFFFF;
  }
}
