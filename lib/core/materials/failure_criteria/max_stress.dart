/// Maximum stress failure criterion for composite plies.
class MaxStressCriterion {
  /// Returns failure mode or null if no failure.
  static FailureMode? failureMode({
    required double sigma1,
    required double sigma2,
    required double tau12,
    required double xt,
    required double xc,
    required double yt,
    required double yc,
    required double s12,
  }) {
    if (sigma1 > 0 && sigma1 >= xt) return FailureMode.longitudinalTension;
    if (sigma1 < 0 && sigma1.abs() >= xc) return FailureMode.longitudinalCompression;
    if (sigma2 > 0 && sigma2 >= yt) return FailureMode.transverseTension;
    if (sigma2 < 0 && sigma2.abs() >= yc) return FailureMode.transverseCompression;
    if (tau12.abs() >= s12) return FailureMode.shear;
    return null;
  }

  static bool fails({
    required double sigma1,
    required double sigma2,
    required double tau12,
    required double xt,
    required double xc,
    required double yt,
    required double yc,
    required double s12,
  }) =>
      failureMode(
        sigma1: sigma1,
        sigma2: sigma2,
        tau12: tau12,
        xt: xt,
        xc: xc,
        yt: yt,
        yc: yc,
        s12: s12,
      ) != null;

  /// Minimum reserve factor across all modes.
  static double reserveFactor({
    required double sigma1,
    required double sigma2,
    required double tau12,
    required double xt,
    required double xc,
    required double yt,
    required double yc,
    required double s12,
  }) {
    final rfList = <double>[];
    if (sigma1 > 0 && xt > 0) rfList.add(xt / sigma1);
    if (sigma1 < 0 && xc > 0) rfList.add(xc / sigma1.abs());
    if (sigma2 > 0 && yt > 0) rfList.add(yt / sigma2);
    if (sigma2 < 0 && yc > 0) rfList.add(yc / sigma2.abs());
    if (tau12.abs() > 0 && s12 > 0) rfList.add(s12 / tau12.abs());
    if (rfList.isEmpty) return double.infinity;
    return rfList.reduce((a, b) => a < b ? a : b);
  }
}

enum FailureMode {
  longitudinalTension,
  longitudinalCompression,
  transverseTension,
  transverseCompression,
  shear,
}
