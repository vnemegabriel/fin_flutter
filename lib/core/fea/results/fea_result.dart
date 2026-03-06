import '../../math/matrix.dart';

/// FEA structural analysis results: natural frequencies and mode shapes.
class FEAResult {
  /// Number of modes extracted.
  final int nModes;

  /// Natural frequencies ω in rad/s.
  final List<double> naturalFrequenciesRadPS;

  /// Natural frequencies in Hz.
  List<double> get naturalFrequenciesHz =>
      naturalFrequenciesRadPS.map((w) => w / (2 * 3.14159265358979)).toList();

  /// Mode shapes: columns of this matrix (n_dof × n_modes).
  final Matrix modeShapes;

  /// Generalized masses: phi_i^T * M * phi_i (should be ~1 if M-normalized).
  final List<double> generalizedMasses;

  const FEAResult({
    required this.nModes,
    required this.naturalFrequenciesRadPS,
    required this.modeShapes,
    required this.generalizedMasses,
  });

  /// Get the i-th mode shape as a list.
  List<double> modeShape(int i) =>
      List.generate(modeShapes.rows, (r) => modeShapes.get(r, i));

  @override
  String toString() =>
      'FEAResult($nModes modes, f1=${naturalFrequenciesHz.first.toStringAsFixed(2)}Hz)';
}
