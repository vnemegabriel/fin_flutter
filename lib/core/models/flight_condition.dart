import 'dart:math';

/// Atmospheric and flight conditions.
class FlightCondition {
  /// Altitude (m).
  final double altitude;

  /// Mach number (dimensionless).
  final double mach;

  /// True airspeed (m/s).
  final double velocity;

  /// Air density (kg/m³).
  final double density;

  /// Dynamic pressure q = 0.5 * rho * V² (Pa).
  double get dynamicPressure => 0.5 * density * velocity * velocity;

  /// Speed of sound (m/s) at given altitude (ISA).
  double get speedOfSound => velocity / (mach > 0 ? mach : 1e-10);

  const FlightCondition({
    required this.altitude,
    required this.mach,
    required this.velocity,
    required this.density,
  });

  /// Create from altitude using ISA standard atmosphere.
  factory FlightCondition.isa({required double altitude, required double mach}) {
    final atm = _ISA.compute(altitude);
    final V = mach * atm.speedOfSound;
    return FlightCondition(
      altitude: altitude,
      mach: mach,
      velocity: V,
      density: atm.density,
    );
  }

  /// Create from velocity and altitude.
  factory FlightCondition.fromVelocity({
    required double altitude,
    required double velocity,
  }) {
    final atm = _ISA.compute(altitude);
    final mach = velocity / atm.speedOfSound;
    return FlightCondition(
      altitude: altitude,
      mach: mach,
      velocity: velocity,
      density: atm.density,
    );
  }

  Map<String, dynamic> toJson() => {
        'altitude': altitude,
        'mach': mach,
        'velocity': velocity,
        'density': density,
      };

  factory FlightCondition.fromJson(Map<String, dynamic> json) => FlightCondition(
        altitude: (json['altitude'] as num).toDouble(),
        mach: (json['mach'] as num).toDouble(),
        velocity: (json['velocity'] as num).toDouble(),
        density: (json['density'] as num).toDouble(),
      );

  @override
  String toString() =>
      'FlightCondition(alt=${altitude.toStringAsFixed(0)}m, '
      'M=${mach.toStringAsFixed(3)}, '
      'V=${velocity.toStringAsFixed(1)}m/s, '
      'rho=${density.toStringAsFixed(4)}kg/m³)';
}

/// ISA 1976 standard atmosphere.
class _ISA {
  final double density;
  final double temperature;
  final double pressure;
  final double speedOfSound;

  const _ISA({
    required this.density,
    required this.temperature,
    required this.pressure,
    required this.speedOfSound,
  });

  static const double _gamma = 1.4;
  static const double _R = 287.058; // J/(kg·K)
  static const double _g0 = 9.80665; // m/s²

  static _ISA compute(double altitude) {
    // Troposphere: 0-11000 m
    // Tropopause: 11000-20000 m
    // Stratosphere-1: 20000-32000 m
    double T, p;
    const T0 = 288.15; // K
    const p0 = 101325.0; // Pa
    const rho0 = 1.225; // kg/m³

    if (altitude <= 11000) {
      const lapseRate = 0.0065; // K/m
      T = T0 - lapseRate * altitude;
      p = p0 * pow(T / T0, _g0 / (lapseRate * _R)).toDouble();
    } else if (altitude <= 20000) {
      const T11 = 216.65;
      const p11 = 22632.1;
      T = T11;
      p = p11 * exp(-_g0 * (altitude - 11000) / (_R * T11));
    } else if (altitude <= 32000) {
      const T20 = 216.65;
      const p20 = 5474.89;
      const lapseRate = -0.001; // K/m (negative = warming)
      T = T20 + lapseRate * (altitude - 20000);
      p = p20 * pow(T / T20, -_g0 / (lapseRate * _R)).toDouble();
    } else {
      // Simplified extension
      const T32 = 228.65;
      const p32 = 868.019;
      const lapseRate = -0.0028;
      T = T32 + lapseRate * (altitude - 32000);
      p = p32 * pow(T / T32, -_g0 / (lapseRate * _R)).toDouble();
    }

    final rho = p / (_R * T);
    final a = sqrt(_gamma * _R * T);
    return _ISA(density: rho, temperature: T, pressure: p, speedOfSound: a);
  }
}
