import 'dart:math';

/// 1D linear interpolation utilities.
class Interpolation {
  /// Linear interpolation given sorted x-values and corresponding y-values.
  static double linear(List<double> xs, List<double> ys, double x) {
    assert(xs.length == ys.length && xs.isNotEmpty);
    if (x <= xs.first) return ys.first;
    if (x >= xs.last) return ys.last;
    // Binary search
    var lo = 0;
    var hi = xs.length - 1;
    while (hi - lo > 1) {
      final mid = (lo + hi) ~/ 2;
      if (xs[mid] <= x) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    final t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
  }

  /// Bilinear interpolation on a regular grid.
  /// [xs] and [ys] define grid axes (sorted).
  /// [zs] is values at grid points: zs[i][j] = z(xs[i], ys[j]).
  static double bilinear(
    List<double> xs,
    List<double> ys,
    List<List<double>> zs,
    double x,
    double y,
  ) {
    // Find x indices
    final ix = _findIndex(xs, x);
    final iy = _findIndex(ys, y);
    final ix1 = min(ix + 1, xs.length - 1);
    final iy1 = min(iy + 1, ys.length - 1);

    final tx = xs[ix1] == xs[ix] ? 0.0 : (x - xs[ix]) / (xs[ix1] - xs[ix]);
    final ty = ys[iy1] == ys[iy] ? 0.0 : (y - ys[iy]) / (ys[iy1] - ys[iy]);

    final z00 = zs[ix][iy];
    final z10 = zs[ix1][iy];
    final z01 = zs[ix][iy1];
    final z11 = zs[ix1][iy1];

    return z00 * (1 - tx) * (1 - ty) +
        z10 * tx * (1 - ty) +
        z01 * (1 - tx) * ty +
        z11 * tx * ty;
  }

  static int _findIndex(List<double> xs, double x) {
    if (x <= xs.first) return 0;
    if (x >= xs.last) return xs.length - 2;
    var lo = 0;
    var hi = xs.length - 1;
    while (hi - lo > 1) {
      final mid = (lo + hi) ~/ 2;
      if (xs[mid] <= x) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }
}
