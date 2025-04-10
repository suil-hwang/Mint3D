#include "MiNTModel_UnitTests.h"
#include "MiNTModel.h"

#include <iostream>

#include "../Optimization/NewtonDescent.h"

#include "MiNTCommonFunc.h"
#include "MiNTEnergy.h"

#include "../mint3D_hook/Serialization.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream> 
#include "../mint3D_hook/date.h"

namespace MiNT3D {

// Flat matrix to vector
Eigen::VectorXd MiNTModel_UnitTests::MatrixToVector(const Eigen::MatrixXd& mat) {
    Eigen::VectorXd vec(mat.rows() * mat.cols());
    for (int i = 0; i < mat.rows(); i++) {
        vec.segment(i * mat.cols(), mat.cols()) = mat.row(i);
    }
    return vec;
}

// Vector to matrix
Eigen::MatrixXd MiNTModel_UnitTests::VectorToMatrix(const Eigen::VectorXd& vec, int cols) {
    Eigen::MatrixXd mat(vec.size() / cols, cols);
    for (int i = 0; i < vec.size() / cols; i++) {
        mat.row(i) = vec.segment(i * cols, cols);
    }
    return mat;
}

// Convert variable to frames, defined on the extended mesh
Eigen::MatrixXd MiNTModel_UnitTests::ConvertVariableToFrames(const Eigen::VectorXd& x) {
    // sanity checks
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    if (!proj_) {
        std::cout << "The projector is not set!" << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    // Eigen::VectorXd full_x = unselection_matrix_ * x + fixed_normal_frames_;

    Eigen::VectorXd full_x;
    proj_->UnprojectVector(x, full_x);
    full_x += fixed_normal_frames_;

    // assert
    assert(full_x.size() % 9 == 0);

    return VectorToMatrix(full_x, 9);
}

// Convert the extended frames to variable
Eigen::VectorXd MiNTModel_UnitTests::ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames) {
    // sanity checks
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return Eigen::VectorXd(0);
    }

    if (!proj_) {
        std::cout << "The projector is not set!" << std::endl;
        return Eigen::VectorXd(0);
    }

    Eigen::VectorXd full_x = MatrixToVector(ext_frames);
    // return selection_matrix_ * full_x;

    Eigen::VectorXd proj_x;
    proj_->ProjectVector(full_x, proj_x);

    return proj_x;
}

// Solve the smooth frames along with boundary normals
bool MiNTModel_UnitTests::SolveTwoTetTest(MiNT3D::SolveParams solve_params, const Eigen::MatrixXd& input_frames,
                                          Eigen::MatrixXd& output_frames) {
    // Step 0: make sure the mesh has been set
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return false;
    }

    // Step 2: set the selection matrix
    std::vector<bool> free_dofs(init_mesh_->nTets() * 9, true);
    fixed_normal_frames_.setZero(9 * init_mesh_->nTets());

    for (int i = 0; i < init_mesh_->nTets(); i++) {
        if (i == 0) {
            for (int j = 0; j < 3; j++) {
                free_dofs[9 * i + j] = false;
            }
            fixed_normal_frames_.segment<3>(9 * i) = input_frames.row(i).segment<3>(0);
        } else {
            fixed_normal_frames_.segment(9 * i, 9) = input_frames.row(i);
            for (int j = 0; j < 9; j++) {
                free_dofs[9 * i + j] = false;
            }
        }
    }

    proj_ = std::make_unique<Projection>(free_dofs);

    CubeCover::TetMeshConnectivity mesh = *init_mesh_;
    Eigen::MatrixXd V = *init_V_;

    // Step 4: Define the objective function and corresponding logging lambda
    MiNTEnergy energy_model(V, mesh, mesh);
    energy_model.Precompute();
    std::cout << energy_model.show_energy_ << std::endl;

    auto log_obj_func = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd cur_ext_frames = ConvertVariableToFrames(x);

        double smoothness_energy =
            energy_model.ComputeSmoothnessEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, false);

        double energy = 0;
        energy += smoothness_energy;

        // std::cout << "mint weight: " << energy_model.weights_.w_mint << std::endl;
        std::cout << "smoothness_energy: " << smoothness_energy << std::endl;

        return energy;
    };

    auto obj_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv = nullptr,
                        Eigen::SparseMatrix<double>* hess = nullptr, bool is_proj = false, double global_scale=0.,
                        bool recompute = false) {
        Eigen::MatrixXd cur_ext_frames = ConvertVariableToFrames(x);
        std::vector<Eigen::Triplet<double>> hess_T;
        Eigen::VectorXd unit_deriv;

        double smoothness_energy =
            energy_model.ComputeSmoothnessEnergyWithTriplets(cur_ext_frames, deriv, hess ? &hess_T : nullptr, is_proj);

        double unit_energy = energy_model.ComputeUnitNormPenaltyWithTriplets(
            cur_ext_frames, deriv ? &unit_deriv : nullptr, hess ? &hess_T : nullptr, is_proj);

        double energy = 0;
        energy += smoothness_energy;
        energy += unit_energy;

        if (deriv) {
            *deriv += unit_deriv;
            Eigen::VectorXd proj_deriv;
            proj_->ProjectVector(*deriv, proj_deriv);
            *deriv = std::move(proj_deriv);
        }

        if (hess) {
            proj_->ProjectMatrix(hess_T);
            int nvars = proj_->ProjDOFs();
            hess->resize(nvars, nvars);
            hess->setFromTriplets(hess_T.begin(), hess_T.end());
        }

        return energy;
    };

    // Step 5: Find the max step size function
    auto find_max_step = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& dir) {
        return 1.0;
    };

    // Step 5: The optimization
    Eigen::VectorXd x0 = ConvertFramesToVariable(input_frames);
    bool to_proj = true;

    // Saninty check gradient and hessian
    OptSolver::TestFuncGradHessian(obj_func, x0);

    // Newton Solver
    bool is_success = OptSolver::NewtonSolver(obj_func, log_obj_func, find_max_step, log_obj_func, x0, solve_params.reg,
                                              solve_params.inner_iter_max_steps, energy_model.weights_.getGlobalScale(),
                                              solve_params.grad_tol, solve_params.xTol, solve_params.fTol, to_proj,
                                              true, solveData.get());
    std::cout << "is success: " << is_success << std::endl;

    // Step 6: Return
    output_frames = ConvertVariableToFrames(x0);

    std::cout << "output frames: " << output_frames.rows() << ", " << output_frames.cols() << std::endl;

    return true;
}

}   // namespace MiNT3D
