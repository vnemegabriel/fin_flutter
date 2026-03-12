#pragma once
/// @file flight_condition.hpp
/// @brief ISA 1976 standard atmosphere and flight condition.
///
/// Reference: docs/theory.md §1 — ISA 1976 Standard Atmosphere.

#include <cmath>
#include <stdexcept>

namespace ff {

/// @brief ISA 1976 atmosphere state at a given altitude.
struct AtmosphereState {
    double temperature;   ///< Static temperature T [K].
    double pressure;      ///< Static pressure p [Pa].
    double density;       ///< Air density ρ [kg/m³].
    double speed_of_sound;///< Speed of sound a [m/s].
};

namespace detail {

/// Physical constants — ISA 1976.
static constexpr double kT0    = 288.15;   ///< Sea-level temperature [K].
static constexpr double kP0    = 101325.0; ///< Sea-level pressure [Pa].
static constexpr double kR     = 287.058;  ///< Specific gas constant [J/(kg·K)].
static constexpr double kg0    = 9.80665;  ///< Gravitational acceleration [m/s²].
static constexpr double kGamma = 1.4;      ///< Heat capacity ratio.
static constexpr double kL     = 0.0065;   ///< Tropospheric lapse rate [K/m].

} // namespace detail

/// @brief Compute ISA 1976 atmosphere state at altitude h [m].
///
/// Layers implemented:
///   0–11 000 m   : troposphere (linear T, power-law p)  — Eq. 1.1–1.2
///   11 000–20 000 m : tropopause (isothermal)           — Eq. 1.3
///   20 000–32 000 m : stratosphere-1 (linear T)         — Eq. 1.4
///   > 32 000 m   : stratosphere-2 (linear T, simplified)
///
/// @param h Geopotential altitude above MSL [m]. Must be ≥ 0.
/// @return AtmosphereState (T, p, ρ, a) at altitude h.
inline AtmosphereState isa1976(double h) {
    using namespace detail;

    double T, p;

    if (h <= 11000.0) {
        // Troposphere — Eq. 1.1
        T = kT0 - kL * h;
        p = kP0 * std::pow(T / kT0, kg0 / (kL * kR));
    } else if (h <= 20000.0) {
        // Tropopause (isothermal) — Eq. 1.3
        constexpr double T11 = 216.65;
        constexpr double p11 = 22632.1;
        T = T11;
        p = p11 * std::exp(-kg0 * (h - 11000.0) / (kR * T11));
    } else if (h <= 32000.0) {
        // Lower stratosphere — Eq. 1.4
        constexpr double T20  = 216.65;
        constexpr double p20  = 5474.89;
        constexpr double L2   = -0.001; // K/m (warming)
        T = T20 + L2 * (h - 20000.0);
        p = p20 * std::pow(T / T20, -kg0 / (L2 * kR));
    } else {
        // Upper stratosphere (simplified extension)
        constexpr double T32  = 228.65;
        constexpr double p32  = 868.019;
        constexpr double L3   = -0.0028;
        T = T32 + L3 * (h - 32000.0);
        p = p32 * std::pow(T / T32, -kg0 / (L3 * kR));
    }

    const double rho = p / (kR * T);
    // Eq. 1.5 — speed of sound
    const double a   = std::sqrt(kGamma * kR * T);

    return {T, p, rho, a};
}

/// @brief Flight condition: altitude, Mach, velocity, density.
struct FlightCondition {
    double altitude; ///< Altitude [m].
    double mach;     ///< Mach number M [-].
    double velocity; ///< True airspeed V [m/s].
    double density;  ///< Air density ρ [kg/m³].

    /// @brief Dynamic pressure q = ½ρV² [Pa].
    double dynamic_pressure() const { return 0.5 * density * velocity * velocity; }

    /// @brief Construct from altitude and Mach number using ISA 1976.
    static FlightCondition from_isa(double altitude, double mach) {
        const auto atm = isa1976(altitude);
        return {altitude, mach, mach * atm.speed_of_sound, atm.density};
    }

    /// @brief Construct from altitude and true airspeed using ISA 1976.
    static FlightCondition from_velocity(double altitude, double velocity) {
        const auto atm = isa1976(altitude);
        return {altitude, velocity / atm.speed_of_sound, velocity, atm.density};
    }
};

} // namespace ff
