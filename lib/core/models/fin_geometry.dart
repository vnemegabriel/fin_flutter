import 'dart:math';

/// Trapezoidal fin geometry definition.
class FinGeometry {
  /// Semi-span (m) - distance from root to tip.
  final double span;

  /// Root chord (m).
  final double rootChord;

  /// Tip chord (m).
  final double tipChord;

  /// Leading-edge sweep length (m) - horizontal offset from root LE to tip LE.
  final double sweepLength;

  /// Fin thickness (m).
  final double thickness;

  /// Number of fins on the rocket (for mass estimation).
  final int finCount;

  const FinGeometry({
    required this.span,
    required this.rootChord,
    required this.tipChord,
    required this.sweepLength,
    required this.thickness,
    this.finCount = 4,
  });

  /// Planform area (m²) - trapezoidal.
  double get planformArea => 0.5 * (rootChord + tipChord) * span;

  /// Aspect ratio AR = b² / S (using full-span: b = 2*span for exposed fin).
  double get aspectRatio => (2.0 * span * span) / planformArea;

  /// Sweep angle at leading edge (radians).
  double get sweepAngleRad => atan2(sweepLength, span);

  /// Sweep angle at leading edge (degrees).
  double get sweepAngleDeg => sweepAngleRad * 180.0 / pi;

  /// Taper ratio λ = c_t / c_r.
  double get taperRatio => tipChord / rootChord;

  /// Mean aerodynamic chord (m).
  double get meanAerodynamicChord =>
      (2.0 / 3.0) *
      rootChord *
      (1.0 + taperRatio + taperRatio * taperRatio) /
      (1.0 + taperRatio);

  /// x-coordinate of leading edge at spanwise location y (0=root, span=tip).
  double leadingEdgeX(double y) => sweepLength * (y / span);

  /// x-coordinate of trailing edge at spanwise location y.
  double trailingEdgeX(double y) => leadingEdgeX(y) + chordAt(y);

  /// Chord length at spanwise station y.
  double chordAt(double y) => rootChord + (tipChord - rootChord) * (y / span);

  Map<String, dynamic> toJson() => {
        'span': span,
        'rootChord': rootChord,
        'tipChord': tipChord,
        'sweepLength': sweepLength,
        'thickness': thickness,
        'finCount': finCount,
      };

  factory FinGeometry.fromJson(Map<String, dynamic> json) => FinGeometry(
        span: (json['span'] as num).toDouble(),
        rootChord: (json['rootChord'] as num).toDouble(),
        tipChord: (json['tipChord'] as num).toDouble(),
        sweepLength: (json['sweepLength'] as num).toDouble(),
        thickness: (json['thickness'] as num).toDouble(),
        finCount: json['finCount'] as int? ?? 4,
      );

  @override
  String toString() =>
      'FinGeometry(span=${span.toStringAsFixed(3)}m, '
      'root=${rootChord.toStringAsFixed(3)}m, '
      'tip=${tipChord.toStringAsFixed(3)}m, '
      'sweep=${sweepAngleDeg.toStringAsFixed(1)}°)';
}
