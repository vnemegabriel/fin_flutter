/// Stability margin and safety factor computation.
class StabilityMargin {
  /// Safety margin factor for flutter (must be > 1).
  static double flutterMargin(double vFlutter, double vMaxFlight) {
    if (vMaxFlight <= 0) return double.infinity;
    return vFlutter / vMaxFlight;
  }

  /// Safety margin factor for divergence.
  static double divergenceMargin(double vDivergence, double vMaxFlight) {
    if (vMaxFlight <= 0) return double.infinity;
    return vDivergence / vMaxFlight;
  }

  /// Returns a safety status based on the flutter margin.
  static SafetyStatus flutterStatus(double margin) {
    if (margin >= 1.5) return SafetyStatus.safe;
    if (margin >= 1.2) return SafetyStatus.warning;
    if (margin >= 1.0) return SafetyStatus.marginal;
    return SafetyStatus.critical;
  }

  /// Recommended minimum flutter margin per MIL-SPEC-9490D: 1.15.
  static const double minimumFlutterMargin = 1.15;

  /// Recommended design flutter margin (20% above flight speed): 1.20.
  static const double designFlutterMargin = 1.20;
}

enum SafetyStatus { safe, warning, marginal, critical }

extension SafetyStatusExt on SafetyStatus {
  String get label => switch (this) {
        SafetyStatus.safe => 'SAFE',
        SafetyStatus.warning => 'WARNING',
        SafetyStatus.marginal => 'MARGINAL',
        SafetyStatus.critical => 'CRITICAL',
      };

  bool get isAcceptable => this == SafetyStatus.safe || this == SafetyStatus.warning;
}
