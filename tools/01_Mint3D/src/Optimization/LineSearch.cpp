
#include "LineSearch.h"
#include <iostream>

namespace OptSolver {
// Backtracking line search with Armijo condition
double BacktrackingArmijo(
    const Eigen::VectorXd& x, const Eigen::VectorXd& grad, const Eigen::VectorXd& dir,
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_func,
    const double global_scale, const double alpha_init) {
    const double c = 0.2;
    const double rho = 0.5;
    double alpha = alpha_init;

    double gradientDirection = grad.dot(dir);
    if (gradientDirection > 0) {
        std::cerr << "Warning: the gradient direction is not a descent direction." << std::endl;
        return alpha = 0.;
    }

    Eigen::VectorXd xNew = x + alpha * dir;
    double fNew = obj_func(xNew, nullptr, nullptr, false, global_scale, false);
    double f = obj_func(x, nullptr, nullptr, false, global_scale, false);
    const double cache = c * grad.dot(dir);

    while (fNew > f + alpha * cache) {
        alpha *= rho;
        xNew = x + alpha * dir;
        fNew = obj_func(xNew, nullptr, nullptr, false, global_scale, false);
    }

    return alpha;
}
}   // namespace OptSolver
