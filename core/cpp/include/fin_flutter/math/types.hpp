#pragma once
/// @file types.hpp
/// @brief Eigen matrix type aliases used throughout fin_flutter.

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ff {

using Matrix3d  = Eigen::Matrix3d;
using MatrixXd  = Eigen::MatrixXd;
using VectorXd  = Eigen::VectorXd;
using Vector3d  = Eigen::Vector3d;
using SparseXd  = Eigen::SparseMatrix<double>;

} // namespace ff
