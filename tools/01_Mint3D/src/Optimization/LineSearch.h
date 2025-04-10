#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace OptSolver {
// Backtracking line search with Armijo condition
double BacktrackingArmijo(
    const Eigen::VectorXd& x, const Eigen::VectorXd& grad, const Eigen::VectorXd& dir,
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_func,
    const double global_scale, const double alpha_init = 1.0);
}   // namespace OptSolver