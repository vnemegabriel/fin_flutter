import 'dart:math';

/// Tsai-Hill failure criterion for composite plies.
class TsaiHillCriterion {
  /// Returns failure index. FI >= 1.0 indicates failure.
  static double failureIndex({
    required double sigma1,
    required double sigma2,
    required double tau12,
    required double xt,
    required double xc,
    required double yt,
    required double yc,
    required double s12,
  }) {
    // Use appropriate strength based on sign of stress
    final X = sigma1 >= 0 ? xt : xc;
    final Y = sigma2 >= 0 ? yt : yc;

    return (sigma1 / X) * (sigma1 / X) -
        (sigma1 * sigma2) / (X * X) +
        (sigma2 / Y) * (sigma2 / Y) +
        (tau12 / s12) * (tau12 / s12);
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
      failureIndex(
        sigma1: sigma1,
        sigma2: sigma2,
        tau12: tau12,
        xt: xt,
        xc: xc,
        yt: yt,
        yc: yc,
        s12: s12,
      ) >= 1.0;

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
    final fi = failureIndex(
      sigma1: sigma1,
      sigma2: sigma2,
      tau12: tau12,
      xt: xt,
      xc: xc,
      yt: yt,
      yc: yc,
      s12: s12,
    );
    if (fi <= 0) return double.infinity;
    return 1.0 / sqrt(fi);
  }
}
