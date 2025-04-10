#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>

#include "MintEnergy.h"

#include "ObjFunc.h"
#include "Projection.h"

// class for objective function
namespace MiNT3D {

// class ObjFunc

// Welcome to the wild and wonderful ObjFuncZoo.h

// This file contains a zoo of objective functions that can be used in optimization problems.

// We generate different objective functiosn and logging functions here.

ObjFunc normalAlignOpt(MiNTEnergy& energy_model, std::unique_ptr<Projection> proj_) {
    auto obj_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv = nullptr,
                        Eigen::SparseMatrix<double>* hess = nullptr, bool is_proj = false, double global_scale) {
        Eigen::VectorXd full_x;
        proj_->UnprojectVector(x, full_x);
        full_x += proj_->fixed_vars_;

        Eigen::MatrixXd cur_ext_frames = proj_->ConvertVariableToFrames(x);
        std::vector<Eigen::Triplet<double>> hess_T;

        std::vector<Eigen::Triplet<double>> unit_hess;
        std::vector<Eigen::Triplet<double>> mint_hess;
        std::vector<Eigen::Triplet<double>> smooth_hess;
        std::vector<Eigen::Triplet<double>> bound_hess;
        std::vector<Eigen::Triplet<double>> sjac_hess;
        std::vector<Eigen::Triplet<double>> det_hess;
        std::vector<Eigen::Triplet<double>> unit_barrier_hess;

        Eigen::VectorXd unit_deriv;
        Eigen::VectorXd integrability_deriv;
        Eigen::VectorXd bound_penalty_deriv;
        // Eigen::VectorXd bound_viscosity_deriv;
        Eigen::VectorXd sjac_deriv;
        Eigen::VectorXd det_deriv;
        Eigen::VectorXd unit_barrier_deriv;

        bool is_krush = energy_model.weights_.b_use_kruskal_tensors_for_sym_smoothness;

        // energy_model.weights_ = energy_model.weights_.setGlobalScale(global_scale);

        double smoothness_energy = energy_model.ComputeSmoothnessEnergyWithTriplets(
            cur_ext_frames, deriv, hess ? &smooth_hess : nullptr, is_proj);
        double mint_energy = energy_model.ComputeIntegrabilityEnergyWithTriplets(
            cur_ext_frames, deriv ? &integrability_deriv : nullptr, hess ? &mint_hess : nullptr, is_proj, is_krush);
        double unit_energy = energy_model.ComputeUnitNormPenaltyWithTriplets(
            cur_ext_frames, deriv ? &unit_deriv : nullptr, hess ? &unit_hess : nullptr, is_proj);
        double boundary_penalty_energy = energy_model.ComputeUnitNormalBoundaryPenaltyWithTriplets(
            cur_ext_frames, deriv ? &bound_penalty_deriv : nullptr, hess ? &bound_hess : nullptr, is_proj);

        // double sjac_energy = energy_model.ComputeScaledJacobianPenaltyWithTriplets(
        //     cur_ext_frames, deriv ? &sjac_deriv : nullptr, hess ? &sjac_hess : nullptr, is_proj);
        double det_energy = energy_model.ComputeDeterminantPenaltyWithTriplets(
            cur_ext_frames, deriv ? &det_deriv : nullptr, hess ? &det_hess : nullptr, is_proj);
        // double unit_barrier_energy = energy_model.ComputeUnitBarrierPenaltyWithTriplets(
        //     cur_ext_frames, deriv ? &unit_barrier_deriv : nullptr, hess ? &unit_barrier_hess : nullptr, is_proj);

        double energy = 0;
        energy += smoothness_energy;
        energy += mint_energy;
        energy += unit_energy;
        energy += boundary_penalty_energy;
        // energy += sjac_energy;
        energy += det_energy;
        // energy += unit_barrier_energy;

        if (deriv) {
            *deriv += integrability_deriv;
            *deriv += unit_deriv;
            *deriv += bound_penalty_deriv;
            // *deriv += sjac_deriv;
            *deriv += det_deriv;
            // *deriv += unit_barrier_deriv;
            Eigen::VectorXd proj_deriv;
            proj_->ProjectVector(*deriv, proj_deriv);
            *deriv = std::move(proj_deriv);
        }

        if (hess) {
            double hsize = unit_hess.size() + mint_hess.size() + smooth_hess.size() + bound_hess.size() +
                           sjac_hess.size() + det_hess.size() + unit_barrier_hess.size();
            hess_T.reserve(hsize);
            hess_T.insert(hess_T.end(), unit_hess.begin(), unit_hess.end());
            hess_T.insert(hess_T.end(), mint_hess.begin(), mint_hess.end());
            hess_T.insert(hess_T.end(), smooth_hess.begin(), smooth_hess.end());
            hess_T.insert(hess_T.end(), bound_hess.begin(), bound_hess.end());
            // hess_T.insert(hess_T.end(), sjac_hess.begin(), sjac_hess.end());
            hess_T.insert(hess_T.end(), det_hess.begin(), det_hess.end());
            // hess_T.insert(hess_T.end(), unit_barrier_hess.begin(), unit_barrier_hess.end());
            proj_->ProjectMatrix(hess_T);
            int nvars = proj_->ProjDOFs();
            hess->resize(nvars, nvars);
            hess->setFromTriplets(hess_T.begin(), hess_T.end());
        }

        return energy;
    };

    ObjFunc obj_func_;
    obj_func_.set_obj_func(obj_func);
    return obj_func_;
}

// ObjFunc createObjFunc(EnergyModel& energyModel, Proj& proj)
// {
//     auto obj_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv = nullptr,
//                         Eigen::SparseMatrix<double>* hess = nullptr, bool is_proj = false, double global_scale) {
//         // The body of the lambda function goes here.
//         // You can use energy_model and proj within this function.
//     };

//     return ObjFunc(obj_func);
// }

}   // namespace MiNT3D
