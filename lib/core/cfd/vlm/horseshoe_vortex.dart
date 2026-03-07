import 'dart:math';

/// 3D vector helper.
class Vec3 {
  final double x, y, z;
  const Vec3(this.x, this.y, this.z);

  Vec3 operator +(Vec3 other) => Vec3(x + other.x, y + other.y, z + other.z);
  Vec3 operator -(Vec3 other) => Vec3(x - other.x, y - other.y, z - other.z);
  Vec3 operator *(double s) => Vec3(x * s, y * s, z * s);

  double dot(Vec3 other) => x * other.x + y * other.y + z * other.z;

  Vec3 cross(Vec3 other) => Vec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x,
      );

  double get norm => sqrt(x * x + y * y + z * z);
  Vec3 get normalized {
    final n = norm;
    if (n < 1e-14) return const Vec3(0, 0, 0);
    return Vec3(x / n, y / n, z / n);
  }

  @override
  String toString() => 'Vec3(${x.toStringAsFixed(4)}, ${y.toStringAsFixed(4)}, ${z.toStringAsFixed(4)})';
}

/// Computes induced velocities using the Biot-Savart law.
class BiotSavart {
  static const double _cutoff = 1e-6; // core radius to avoid singularity

  /// Velocity induced by a finite vortex segment A→B with circulation Gamma.
  static Vec3 finiteSegment({
    required Vec3 P,   // field point
    required Vec3 A,   // segment start
    required Vec3 B,   // segment end
    required double gamma, // circulation strength (m²/s)
  }) {
    final r1 = P - A;
    final r2 = P - B;
    final r1xr2 = r1.cross(r2);
    final crossMagSq = r1xr2.dot(r1xr2);

    if (crossMagSq < _cutoff * _cutoff) {
      return const Vec3(0, 0, 0); // point on or near vortex axis
    }

    final r1Mag = r1.norm;
    final r2Mag = r2.norm;

    if (r1Mag < _cutoff || r2Mag < _cutoff) {
      return const Vec3(0, 0, 0);
    }

    final factor = gamma / (4.0 * pi) *
        (1.0 / r1Mag + 1.0 / r2Mag) *
        (r1.dot(r2) / (r1Mag * r2Mag) - 1.0);
    // Wait, let me use the standard formula:
    // V = (Gamma/4pi) * (r1 x r2) / |r1 x r2|^2 * ((r1/|r1| - r2/|r2|) . (r1-r2)/|r1-r2|)
    // Simplified form:
    // V = Gamma/(4pi) * (r1 x r2) / |r1 x r2|^2 * (1/r1 + 1/r2) * cos(theta1+theta2)/2 ...
    // Standard correct form:
    final AB = B - A;
    final cosTerm = AB.normalized.dot(r1.normalized) - AB.normalized.dot(r2.normalized);
    final strength = gamma / (4.0 * pi) * cosTerm / crossMagSq;

    return r1xr2 * strength;
  }

  /// Velocity induced by a semi-infinite vortex from point A in direction +x (wake).
  static Vec3 semiInfinite({
    required Vec3 P,
    required Vec3 A,
    required double gamma,
    required Vec3 direction, // unit vector of wake direction
  }) {
    final r = P - A;
    final rCrossd = r.cross(direction);
    final crossMagSq = rCrossd.dot(rCrossd);

    if (crossMagSq < _cutoff * _cutoff) {
      return const Vec3(0, 0, 0);
    }

    final rDotd = r.dot(direction);
    final rMag = r.norm;
    final factor = gamma / (4.0 * pi) * (1.0 + rDotd / rMag) / crossMagSq;

    return rCrossd * factor;
  }
}
