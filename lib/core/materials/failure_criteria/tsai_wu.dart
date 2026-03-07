/// Tsai-Wu failure criterion for composite plies.
class TsaiWuCriterion {
  /// Evaluate Tsai-Wu failure index for given stress state.
  ///
  /// Returns the failure index FI. FI >= 1.0 indicates failure.
  /// Reserve factor RF = 1 / FI.
  ///
  /// [sigma1] - longitudinal stress (Pa)
  /// [sigma2] - transverse stress (Pa)
  /// [tau12]  - in-plane shear stress (Pa)
  /// Strength values all in Pa (positive = tensile).
  static double failureIndex({
    required double sigma1,
    required double sigma2,
    required double tau12,
    required double xt, // longitudinal tensile strength
    required double xc, // longitudinal compressive strength (positive)
    required double yt, // transverse tensile strength
    required double yc, // transverse compressive strength (positive)
    required double s12, // shear strength
  }) {
    final f1 = 1.0 / xt - 1.0 / xc;
    final f2 = 1.0 / yt - 1.0 / yc;
    final f11 = 1.0 / (xt * xc);
    final f22 = 1.0 / (yt * yc);
    final f66 = 1.0 / (s12 * s12);
    final f12 = -0.5 * (f11 * f22).abs() == 0
        ? 0.0
        : -0.5 * (f11 * f22 < 0 ? -((-f11 * f22).abs()) : (f11 * f22).abs());

    return f1 * sigma1 +
        f2 * sigma2 +
        f11 * sigma1 * sigma1 +
        f22 * sigma2 * sigma2 +
        f66 * tau12 * tau12 +
        2.0 * f12 * sigma1 * sigma2;
  }

  /// Returns true if the ply has failed.
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

  /// Reserve factor (RF = 1/FI; RF >= 1 means safe).
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
    return 1.0 / fi;
  }
}
