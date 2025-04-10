#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

namespace OptSolver {

// Inherit from this and implement functions and store data for a specific problem.
// class EnergyData
// {
// public:
// 	EnergyData() {}
// 	EnergyData(const Eigen::VectorXd& x, double energy, const Eigen::VectorXd& grad, const
// Eigen::SparseMatrix<double>& hessian, bool is_proj) 		: x_(x), energy_(energy), grad_(grad),
// hessian_(hessian), is_proj_(is_proj) {}

// 	// virtual void log_energy_parts_at(Eigen::VectorXd& x0) = 0;

// 	Eigen::VectorXd x_;
// 	double energy_;
// 	Eigen::VectorXd grad_;
// 	Eigen::SparseMatrix<double> hessian_;
// 	bool is_proj_;

// };

class SolveData {
   public:
    std::vector<double> total_time;
    std::vector<double> assembly_time;
    std::vector<double> solve_time;
    std::vector<double> line_search_time;

    std::vector<double> identity_weight;

    std::vector<double> total_energy;
    std::vector<double> energy_diff;
    std::vector<double> solve_residual;
    std::vector<double> gradient_norm;
    std::vector<double> gradient_norm_step_start;

    std::vector<double> global_scale_log;
    

    bool logToFile(std::string logPath);
};

// Newton solver with line search
/*
 * Input:
 * obj_func: the objective function, which takes x as input and returns the function value, gradient, and hessian,
 * together with a boolean indicating whether the hessian matrix is PSD projected find_max_step: the function to find
 * the maximum step size, which takes x, direction, and returns the maximum step size x0: the initial guess num_iter:
 * the maximum number of iterations grad_tol: the termination tolerance of the gradient x_tol: the tolerance of the
 * solution f_tol: the tolerance of the function value display_info: whether to display the information
 *
 */
bool NewtonSolver(
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_func,
    std::function<double(const Eigen::VectorXd&)> log_obj_func,
    std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> find_max_step,
    std::function<double(const Eigen::VectorXd&)> per_inner_iter_update, Eigen::VectorXd& x0,
    double& reg, int num_iter = 1000, double global_scale = 1., double grad_tol = 1e-14, double x_tol = 0,
    double f_tol = 0, bool is_proj = false, bool display_info = false, bool log_inner_state = false,
    SolveData* solve_data = nullptr);

// Alternative linear solvers
int CallNASOQ(Eigen::SparseMatrix<double, Eigen::ColMajor, int> A, Eigen::Matrix<double, Eigen::Dynamic, 1> b,
              Eigen::Matrix<double, Eigen::Dynamic, 1>& x);
// BUGGED
void CallSPQR(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& neg_grad, Eigen::VectorXd& delta_x);
// DOESN'T ACTUALLY HELP CONVERGENCE
void CallGMRES(const Eigen::SparseMatrix<double>& H, const Eigen::VectorXd& neg_grad, Eigen::VectorXd& delta_x);

// Test the function gradient and hessian
void TestFuncGradHessian(
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_func,
    const Eigen::VectorXd& x0);

}   // namespace OptSolver