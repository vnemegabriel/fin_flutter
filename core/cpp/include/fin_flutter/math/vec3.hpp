#pragma once
/// @file vec3.hpp
/// @brief Lightweight 3-D vector for Biot-Savart and panel geometry.
///
/// Intentionally separate from Eigen::Vector3d to keep the CFD
/// panel-building code readable without Eigen boilerplate.

#include <cmath>
#include <string>

namespace ff {

/// @brief 3-D Cartesian vector (double precision).
struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)      const { return {x*s,   y*s,   z*s};   }

    /// @brief Dot product.
    double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }

    /// @brief Cross product.
    Vec3 cross(const Vec3& o) const {
        return {y*o.z - z*o.y,
                z*o.x - x*o.z,
                x*o.y - y*o.x};
    }

    /// @brief Euclidean norm.
    double norm() const { return std::sqrt(x*x + y*y + z*z); }

    /// @brief Unit vector; returns zero vector if norm < 1e-14.
    Vec3 normalized() const {
        const double n = norm();
        return (n < 1e-14) ? Vec3{} : Vec3{x/n, y/n, z/n};
    }
};

} // namespace ff
