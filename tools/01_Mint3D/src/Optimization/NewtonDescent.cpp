#include <fstream>
#include <iomanip>
#include <iostream>

#include <Eigen/CholmodSupport>
// #include <Eigen/SPQRSupport>
#include <Eigen/Sparse>

#include <cholmod.h>   // Include CHOLMOD header
#include <SuiteSparseQR.hpp>

#include <Eigen/IterativeLinearSolvers>
// #include <Eigen/Unsupported>
#include <unsupported/Eigen/IterativeSolvers>
#include "CholeskyPreconditioner.h"
#include "LineSearch.h"
#include "NewtonDescent.h"

#include <Eigen/Sparse>

#include "../mint3D_hook/Serialization.h"

// #include <nasoq_eigen.h>
// #include <nasoq/nasoq.h>
// #include <nasoq/nasoq_eigen.h>

#include "Timer.h"

namespace OptSolver {
// Newton solver with line search

bool NewtonSolver(
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_func,
    std::function<double(const Eigen::VectorXd&)> log_obj_func,
    std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> find_max_step,
    std::function<double(const Eigen::VectorXd&)> per_inner_iter_update, Eigen::VectorXd& x0, double& reg, int num_iter,
    double global_scale, double grad_tol, double xTol, double fTol, bool is_proj, bool display_info,
    bool log_inner_state, SolveData* solve_data) {
    const int DIM = x0.rows();
    // Eigen::VectorXd randomVec = x0;
    // randomVec.setRandom();
    // x0 += 1e-6 * randomVec;
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(DIM);
    Eigen::VectorXd grad_step_start = Eigen::VectorXd::Zero(DIM);

    Eigen::SparseMatrix<double> hessian;

    Eigen::VectorXd neg_grad, delta_x;
    double max_step_size = 1.0;
    double reg_inp = reg;
    // double reg = 1e-8;

    // is_proj = !is_proj;
    Timer<std::chrono::high_resolution_clock> total_timer;
    double total_assembling_time = 0;
    double total_solving_time = 0;
    double total_linesearch_time = 0;

    total_timer.start();
    std::ofstream optInfo;
    if (display_info) {
        std::cout << "gradient tol: " << grad_tol << ", function update tol: " << fTol
                  << ", variable update tol: " << xTol << ", maximum iteration: " << num_iter << std::endl
                  << std::endl;
    }
    int i = 0;

    double f = obj_func(x0, nullptr, nullptr, false, 1., false);
    if (f == 0) {
        std::cout << "energy = 0, return" << std::endl;
        return true;
    }

    Eigen::SparseMatrix<double> I(DIM, DIM);
    I.setIdentity();

    bool is_small_perturb_needed = true;

    double max_grad_prev = 1;
    double max_grad_curr = 1;

    double alpha_line_search = 1.0;

    int toggle = 0;
    int toggle_max = 6;

    double curr_grad_tol = std::max(grad_tol * global_scale, 1e-12);
    double func_progress = 0;

    double proj_reg = 1e-12;

    per_inner_iter_update(x0);

    try {
        for (; i < num_iter; i++) {
            if (display_info) {
                // std::cout << "\niter: " << i << std::endl;
                std::cout << "\nInner Step: " << i << std::endl;
            }

            if (log_inner_state) {
                log_obj_func(x0);
            }
            Eigen::VectorXd xprev = x0;

            double current_scale = 1.;   // this is bugged if not set to 1!!! TODO: probably just remove this...

            Timer<std::chrono::high_resolution_clock> local_timer;
            local_timer.start();
            double f = obj_func(x0, &grad_step_start, &hessian, is_proj, current_scale, true);
            grad = grad_step_start;
            local_timer.stop();
            double localAssTime = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
            total_assembling_time += localAssTime;

            local_timer.start();

            Eigen::SparseMatrix<double> H = hessian;
            std::cout << "num of nonzeros: " << H.nonZeros() << ", rows: " << H.rows() << ", cols: " << H.cols()
                      << ", Sparsity: " << H.nonZeros() * 100.0 / (H.rows() * H.cols()) << "%, ||H||: " << H.norm()
                      << std::endl;

            if (is_small_perturb_needed && is_proj) {   // is_small_perturb_needed &&
                // due to the numerical issue, we may need to add a small perturbation to the PSD projected hessian
                // matrix
                H += proj_reg * I;
                std::cout << "using proj hessian with perturb" << std::endl;
            }

            if (is_proj) {
                std::cout << "using proj hessian" << std::endl;
            } else {
                H += reg * I;
                std::cout << "using diagonally regualized hessian" << std::endl;
            }

            Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver(H);

            //        bool increase_reg = false;
            bool switched = false;
            double rate = 0.;
            double prev_reg = reg;

            auto inf_norm = [](const Eigen::SparseMatrix<double>& mat) {
                double lInfNorm = 0;
                for (int k = 0; k < mat.outerSize(); ++k) {
                    double rowSum = 0;
                    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                        rowSum += abs(it.value());
                    }
                    lInfNorm = std::max(lInfNorm, rowSum);
                }
                return lInfNorm;
            };

            while (solver.info() != Eigen::Success || rate == 0) {
                bool use_qr_solver = false;
                if (solver.info() != Eigen::Success) {
                    if (display_info) {
                        if (is_proj) {
                            std::cout << "some small perturb is needed to remove round-off error, current reg = " << reg
                                      << ", ||H|| = " << H.norm() << std::endl;
                            is_proj = false;
                        } else {
                            std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
                            if (reg > 1e30) {
                                return false;
                            }
                        }
                    }

                    if (is_proj) {
                        is_small_perturb_needed = true;
                    }

                    //                increase_reg = true;
                    reg = std::max(2 * reg, 1e-17);

                    if (is_proj) {
                        H = H + proj_reg * I;
                    } else {
                        H = hessian + reg * I;
                    }

                    solver.compute(H);

                    if (reg > 1e8 && !is_proj) {
                        std::cout << "reg is too large, try to switch to use SPD hessian instead." << std::endl;

                        Eigen::SparseMatrix<double> tmp_hess;
                        double f = obj_func(x0, &grad, &tmp_hess, true, current_scale, false);

                        if (inf_norm(tmp_hess - hessian) > reg) {
                            std::cout << "PSD hessian are further away from the actual one, use diagonal one"
                                      << std::endl;
                        } else {
                            std::cout << "switch to use PSD hessian" << std::endl;
                            is_proj = true;
                            hessian = std::move(tmp_hess);
                            H = hessian + reg * I;   // reg ??? proj_reg
                            solver.compute(H);
                        }
                    }

                    std::cout << "make it here " << std::endl;
                }
                if (solver.info() == Eigen::Success) {
                    // use qr if problem has large residual

                    neg_grad = -grad;
                    delta_x = solver.solve(neg_grad);

                    double residual = (H * delta_x + grad).norm();

                    local_timer.stop();
                    double local_solving_time = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
                    total_solving_time += local_solving_time;

                    max_step_size = find_max_step(x0, delta_x);

                    local_timer.start();
                    rate = BacktrackingArmijo(x0, grad, delta_x, obj_func, current_scale, max_step_size);
                    local_timer.stop();
                    double local_linesearch_time = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
                    total_linesearch_time += local_linesearch_time;
                    alpha_line_search = rate;

                    if (rate == 0) {
                        std::cout << "not a descent direction, switch proj state" << std::endl;
                        is_proj = !is_proj;
                        switched = true;
                    }
                    // std::cout << "toggle " << toggle << std::endl;

                    // Update kinetic term and exponential term.
                    if (!is_proj) {
                        std::cout << "adjust reg " << reg << std::endl;


                        // Do line search
                        reg *= 0.5;
                        reg = std::max(reg, 1e-16);
                        toggle = 0;
                        std::cout << "decreased " << reg << std::endl;

                    } else {
                        // Exit condition for projection state.
                        double fnew = obj_func(x0 + rate * delta_x, nullptr, nullptr, is_proj, current_scale, false);
                        //                    if ((((rate < .4) && grad.norm() < .005) || (rate > .999 && grad.norm() <
                        //                    .001) ||
                        //                         ((f - fnew < 1e-3))) &&
                        //                        !switched)   // delta_x.norm() < 1e-8) )

                        bool not_at_start = solve_data->total_time.size() > 10;
                        if ((grad.norm() < 1e-3 || (f - fnew) / fnew < 1e-5) && !switched && not_at_start ||
                            rate < 1e-5 && not_at_start) {
                            is_proj = false;
                            // reg = 1e-12;
                            toggle = 0;

                            double gradnorm = grad.norm();
                            double fprog = (f - fnew) / fnew;
                            std::cout << "Switch back to diag hessian: diag reg " << reg << " f prog " << fprog
                                      << " grad norm " << gradnorm << std::endl;
                        }
                        // reg = 1e-8;
                    }

                    // switched = false;
                }
            }

            double residual = (H * delta_x + grad).norm();

            // Log in the scale where smoothness = 1.
            f = obj_func(x0, &grad, nullptr, is_proj, current_scale, false);

            x0 = x0 + rate * delta_x;

            double fnew = obj_func(x0, &grad, nullptr, is_proj, current_scale, false);
            double func_progress = (f - fnew);

            per_inner_iter_update(xprev);

            if (display_info) {
                std::cout << "line search rate : " << rate << ", actual hessian : " << !is_proj << ", reg = " << reg
                          << std::endl;
                std::cout << "f_old: " << f << ", f_new: " << fnew << ", grad step start: " << grad_step_start.norm()
                          << ", grad step end: " << grad.norm() << ", delta x: " << rate * delta_x.norm()
                          << ", delta_f: " << f - fnew << ", cur global scale: " << global_scale
                          << ", residual: " << residual << std::endl;
                std::cout << "timing info (in total seconds): " << std::endl;
                std::cout << "assembling took: " << total_assembling_time << ", LLT solver took: " << total_solving_time
                          << ", line search took: " << total_linesearch_time << std::endl;
                // std::cout << "reg: " << reg << ", kinetic term: " << kinetic_term << ", total reg: " << reg +
                // kinetic_term_weight * kinetic_term << std::endl;
            }

            if (solve_data) {
                solve_data->total_time.push_back(total_timer.elapsed<std::chrono::milliseconds>() * 1e-3);
                solve_data->assembly_time.push_back(total_assembling_time);
                solve_data->solve_time.push_back(total_solving_time);
                solve_data->line_search_time.push_back(total_linesearch_time);
                solve_data->identity_weight.push_back(reg);
                solve_data->total_energy.push_back(fnew);
                solve_data->energy_diff.push_back(f - fnew);
                solve_data->solve_residual.push_back(residual);
                solve_data->gradient_norm.push_back(grad.norm());
                solve_data->gradient_norm_step_start.push_back(grad_step_start.norm());
                solve_data->global_scale_log.push_back(global_scale);
            }

            if (grad.norm() < curr_grad_tol) {   //  && rate > .999

                if (grad.norm() < 1e-14 && fnew > 1e8) {
                    x0 += 1e-6 * Eigen::VectorXd::Random(x0.rows());
                } else {
                    std::cout << "terminate with gradient L2-norm = " << grad.norm() << std::endl;
                    break;
                }
            }

            // if (!is_proj && prev_reg == reg && i > 25 &&
            //     (rate < std::max(current_scale, 1e-4) || reg < 1e-15)) {   //  && rate > .999
            //     std::cout << "exit inner iteration because making slow progress " << grad.norm() << std::endl;
            //     break;
            // }

            // if (delta_x.norm() < xTol * rate) {
            //     std::cout << "terminate with small variable change, gradient L2-norm = " << grad.norm() << std::endl;
            //     break;
            // }

            if (f - fnew < fTol && i > 1) {
                std::cout << "terminate with small energy change, gradient L2-norm = " << grad.norm() << std::endl;

                if (reg > 1e-16)
                {
                    // reg *= 1e-4;
                    reg = 1e-16;
                }
                break;
            }
        }

    }   // end of try
    catch (const std::runtime_error& e) {
        // Handle the specific exception
        if (std::string(e.what()) == "Energy is nan") {
            std::cerr << "Caught exception: " << e.what() << '\n';
            // Add any additional handling code here
        } else {
            throw;   // Re-throw the exception if it's not the one we're looking for
        }
    } catch (...) {
        // Handle all other exceptions
        throw;
    }

    if (i >= num_iter) {
        std::cout << "terminate with reaching the maximum iteration, with gradient L2-norm = " << grad.norm()
                  << std::endl;
    }

    f = obj_func(x0, &grad, nullptr, false, 1., false);
    std::cout << "end up with energy: " << f << ", gradient: " << grad.norm() << std::endl;

    total_timer.stop();
    if (display_info) {
        std::cout << "total time cost (s): " << total_timer.elapsed<std::chrono::milliseconds>() * 1e-3
                  << ", within that, assembling took: " << total_assembling_time
                  << ", LLT solver took: " << total_solving_time << ", line search took: " << total_linesearch_time
                  << std::endl;
    }

    return true;
}


// Test the function gradient and hessian
void TestFuncGradHessian(
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double, bool)>
        obj_Func,
    const Eigen::VectorXd& x0) {
    Eigen::VectorXd dir = x0;
    dir(0) = 0;
    dir.setRandom();

    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> H;

    double f = obj_Func(x0, &grad, &H, false, 1., false);
    std::cout << "f: " << f << std::endl;
    if (f == 0) return;

    std::cout << "\nFull Objective Function Finite Difference Test" << std::endl;

    for (int i = 3; i < 10; i++) {
        double eps = std::pow(0.1, i);
        Eigen::VectorXd x = x0 + eps * dir;
        Eigen::VectorXd grad1;
        double f1 = obj_Func(x, &grad1, nullptr, false, 1., false);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "energy-gradient: " << (f1 - f) / eps - grad.dot(dir) << std::endl;
        std::cout << "gradient-hessian: " << ((grad1 - grad) / eps - H * dir).norm() << std::endl;
    }
}

bool SolveData::logToFile(std::string logPath) {
    Serialization::writeCSV(total_time, logPath + "/total_time.csv");
    Serialization::writeCSV(assembly_time, logPath + "/assembly_time.csv");
    Serialization::writeCSV(solve_time, logPath + "/solve_time.csv");
    Serialization::writeCSV(line_search_time, logPath + "/line_search_time.csv");

    Serialization::writeCSV(identity_weight, logPath + "/identity_weight.csv");

    Serialization::writeCSV(total_energy, logPath + "/total_energy.csv");
    Serialization::writeCSV(energy_diff, logPath + "/energy_diff.csv");
    Serialization::writeCSV(solve_residual, logPath + "/solve_residual.csv");
    Serialization::writeCSV(gradient_norm, logPath + "/gradient_norm.csv");
    Serialization::writeCSV(gradient_norm_step_start, logPath + "/gradient_norm_step_start.csv");

    Serialization::writeCSV(global_scale_log, logPath + "/global_scale_log.csv");

    return true;
}

}   // namespace OptSolver
