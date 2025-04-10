#include "MiNTModel.h"

#include <iostream>

#include "../Optimization/NewtonDescent.h"

#include "MiNTCommonFunc.h"
#include "MiNTEnergy.h"

#include "../mint3D_hook/Serialization.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace MiNT3D {

// Fix this add in the missing terms.
bool EnergyLog::logToFile(std::string logPath) {
    Serialization::writeCSV(smoothness, logPath + "/smoothness.csv");
    Serialization::writeCSV(unit_penalty, logPath + "/unit_penalty.csv");
    Serialization::writeCSV(symmetric_integrability, logPath + "/symmetric_integrability.csv");
    Serialization::writeCSV(primal_integrability, logPath + "/primal_integrability.csv");
    Serialization::writeCSV(scaled_jacobian, logPath + "/scaled_jacobian.csv");
    Serialization::writeCSV(unit_barrier, logPath + "/unit_barrier.csv");
    Serialization::writeCSV(curl_weight, logPath + "/curl_weight.csv");
    Serialization::writeCSV(asap_combed_smoothness, logPath + "/asap_combed_smoothness.csv");
    Serialization::writeCSV(aiap_combed_smoothness, logPath + "/aiap_combed_smoothness.csv");

    Serialization::writeCSV(combed_smoothness, logPath + "/combed_smoothness.csv");
    Serialization::writeCSV(combed_integrability, logPath + "/combed_integrability.csv");
    Serialization::writeCSV(total_energy_unscale, logPath + "/total_energy_unscale.csv");

    Serialization::writeCSV(global_scale_outiter, logPath + "/global_scale_outer_iter.csv");

    return true;
}

Projection::Projection(const std::vector<bool>& keep_dofs) {
    int fulldofs = keep_dofs.size();
    dofmap.resize(fulldofs);
    constraint_mask_ = Eigen::VectorXd::Zero(fulldofs);
    int idx = 0;
    for (int i = 0; i < fulldofs; i++) {
        if (keep_dofs[i]) {
            dofmap[i] = idx;
            invdofmap.push_back(i);
            idx++;
        } else {
            dofmap[i] = -1;
            constraint_mask_[i] = 1;
        }
    }
}

void Projection::setFixedVars(const Eigen::MatrixXd& ext_frames) {
    // multiply coefficient wise with the constraint mask

    fixed_vars_ = MiNTModel::MatrixToVector(ext_frames);

    fixed_vars_ = fixed_vars_.cwiseProduct(constraint_mask_);
}

void Projection::ProjectVector(const Eigen::VectorXd& full_vec, Eigen::VectorXd& proj_vec) const {
    int projdofs = invdofmap.size();
    proj_vec.resize(projdofs);

    for (int i = 0; i < projdofs; i++) {
        proj_vec[i] = full_vec[invdofmap[i]];
    }
}

void Projection::UnprojectVector(const Eigen::VectorXd& proj_vec, Eigen::VectorXd& full_vec) const {
    int fulldofs = dofmap.size();
    full_vec.resize(fulldofs);
    for (int i = 0; i < fulldofs; i++) {
        full_vec[i] = (dofmap[i] == -1 ? 0.0 : proj_vec[dofmap[i]]);
    }
}

void Projection::ProjectMatrix(std::vector<Eigen::Triplet<double>>& mat) const {
    int dim = mat.size();
    for (int i = 0; i < dim; i++) {
        int r = mat[i].row();
        int c = mat[i].col();
        int pr = dofmap[r];
        int pc = dofmap[c];
        if (pr != -1 && pc != -1) {
            mat[i] = {pr, pc, mat[i].value()};
        } else {
            mat[i] = {0, 0, 0.0};
        }
    }
}

// Add additional layers for the boundary
// The extended mesh is the mesh with additional layers for the boundary by flipping to properly handle the boundary
// conditions
void MiNTModel::AddAdditionalLayerForBoundary() {
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return;
    }

    // Add additional layers for the boundary
    std::vector<int> bnd_faces;
    for (int i = 0; i < init_mesh_->nFaces(); i++) {
        if (init_mesh_->isBoundaryFace(i)) {
            bnd_faces.push_back(i);
        }
    }

    extended_V_.resize(init_V_->rows() + bnd_faces.size(), init_V_->cols());
    extended_V_.block(0, 0, init_V_->rows(), init_V_->cols()) = *init_V_;

    Eigen::MatrixXi extended_T(init_mesh_->nTets() + bnd_faces.size(), 4);
    for (int i = 0; i < init_mesh_->nTets(); i++) {
        extended_T.row(i) << init_mesh_->tetVertex(i, 0), init_mesh_->tetVertex(i, 1), init_mesh_->tetVertex(i, 2),
            init_mesh_->tetVertex(i, 3);
    }

    int n_ext_verts = init_V_->rows();
    int n_ext_tets = init_mesh_->nTets();

    fixed_normal_frames_.setZero(9 * (extended_T.rows()));
    fixed_interior_frames_.setZero(9 * (extended_T.rows()));

    for (int i = 0; i < bnd_faces.size(); i++) {
        int fid = bnd_faces[i];
        int tet_id = init_mesh_->faceTet(fid, 0) != -1 ? init_mesh_->faceTet(fid, 0) : init_mesh_->faceTet(fid, 1);
        // compute the normal of the boundary face
        Eigen::Vector3d v0 = init_V_->row(init_mesh_->faceVertex(fid, 0));
        Eigen::Vector3d v1 = init_V_->row(init_mesh_->faceVertex(fid, 1));
        Eigen::Vector3d v2 = init_V_->row(init_mesh_->faceVertex(fid, 2));

        // project the face opposite vertex to the normal direction
        int opp_vid = -1;
        for (int j = 0; j < 4; j++) {
            int vid = init_mesh_->tetVertex(tet_id, j);
            if (vid != init_mesh_->faceVertex(fid, 0) && vid != init_mesh_->faceVertex(fid, 1) &&
                vid != init_mesh_->faceVertex(fid, 2)) {
                opp_vid = vid;
                break;
            }
        }
        assert(opp_vid != -1);
        Eigen::Vector3d p = init_V_->row(opp_vid);

        // compute the projected point: p0 which can be expressed as p0 = v0 + a1 * (v1 - v0) + a2 * (v2 - v0)
        // and satisfies: (p - p0).dot(v1 - v0) = 0 and (p - p0).dot(v2 - v0) = 0
        Eigen::Matrix2d mat;
        mat << (v1 - v0).squaredNorm(), (v1 - v0).dot(v2 - v0), (v1 - v0).dot(v2 - v0), (v2 - v0).squaredNorm();
        Eigen::Vector2d rhs((p - v0).dot(v1 - v0), (p - v0).dot(v2 - v0));
        Eigen::Vector2d sol = mat.inverse() * rhs;

        Eigen::Vector3d p0 = v0 + sol[0] * (v1 - v0) + sol[1] * (v2 - v0);

        // Eigen::Vector3d flipped_p = 2 * p0 - p;
        Eigen::Vector3d flipped_p =
            (p0 - p) + (v0 + v1 + v2) / 3;   // make new flipped point be relative to face centroid

        extended_V_.row(n_ext_verts) = flipped_p.transpose();

        extended_T.row(n_ext_tets) << init_mesh_->faceVertex(fid, 0), init_mesh_->faceVertex(fid, 1),
            init_mesh_->faceVertex(fid, 2), n_ext_verts;

        // verify it oriented correctly
        Eigen::Matrix3d A;
        A.col(0) = v0 - flipped_p;
        A.col(1) = v1 - flipped_p;
        A.col(2) = v2 - flipped_p;

        Eigen::Vector3d tv0 = init_V_->row(init_mesh_->tetVertex(tet_id, 0));
        Eigen::Vector3d tv1 = init_V_->row(init_mesh_->tetVertex(tet_id, 1));
        Eigen::Vector3d tv2 = init_V_->row(init_mesh_->tetVertex(tet_id, 2));
        Eigen::Vector3d tv3 = init_V_->row(init_mesh_->tetVertex(tet_id, 3));

        Eigen::Matrix3d tetVol;
        tetVol.col(0) = tv1 - tv0;
        tetVol.col(1) = tv2 - tv0;
        tetVol.col(2) = tv3 - tv0;

        double tetVolSign = tetVol.determinant();

        if (A.determinant() * tetVolSign > 0) {
            extended_T.row(n_ext_tets) << init_mesh_->faceVertex(fid, 0), init_mesh_->faceVertex(fid, 1), n_ext_verts,
                init_mesh_->faceVertex(fid, 2);
        }

        // fixed normals
        Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0);
        normal.normalize();

        fixed_normal_frames_.segment<3>(9 * n_ext_tets) = normal;

        n_ext_verts++;
        n_ext_tets++;
    }
    extended_mesh_ = CubeCover::TetMeshConnectivity(extended_T);

    if (!extended_mesh_.isFaceConnected()) {
        throw std::runtime_error("Input mesh not valid.  Probably has inverted tets.");
    }

    extended_T_ = extended_T;

    // for(int i = 0; i < extended_mesh_.nFaces(); i++)
    // {
    //     std::cout << "Face " << i << " is boundary: " << extended_mesh_.isBoundaryFace(i) << std::endl;
    // }
}

// Set the selection matrix
Eigen::SparseMatrix<double> MiNTModel::ComputeSelectionMatrix(const std::vector<bool>& fixed_dofs) {
    std::vector<Eigen::Triplet<double>> T;
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return Eigen::SparseMatrix<double>(0, 0);
    }

    int full_dofs = fixed_dofs.size();
    int rows = 0;
    for (int i = 0; i < full_dofs; i++) {
        // free variable
        if (!fixed_dofs[i]) {
            T.push_back(Eigen::Triplet<double>(rows, i, 1));
            rows++;
        }
    }
    Eigen::SparseMatrix<double> select_mat(rows, full_dofs);
    select_mat.setFromTriplets(T.begin(), T.end());
    return select_mat;
}

// Flat matrix to vector
Eigen::VectorXd MiNTModel::MatrixToVector(const Eigen::MatrixXd& mat) {
    Eigen::VectorXd vec(mat.rows() * mat.cols());
    for (int i = 0; i < mat.rows(); i++) {
        vec.segment(i * mat.cols(), mat.cols()) = mat.row(i);
    }
    return vec;
}

// Vector to matrix
Eigen::MatrixXd MiNTModel::VectorToMatrix(const Eigen::VectorXd& vec, int cols) {
    Eigen::MatrixXd mat(vec.size() / cols, cols);
    for (int i = 0; i < vec.size() / cols; i++) {
        mat.row(i) = vec.segment(i * cols, cols);
    }
    return mat;
}

// Convert variable to frames, defined on the extended mesh
Eigen::MatrixXd MiNTModel::ConvertVariableToFrames(const Eigen::VectorXd& x) {
    // sanity checks
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    if (extended_V_.rows() == 0) {
        std::cerr << "The extended mesh is not set!" << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    if (!proj_) {
        std::cout << "The projector is not set!" << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    // Eigen::VectorXd full_x = unselection_matrix_ * x + fixed_normal_frames_;

    Eigen::VectorXd full_x;
    proj_->UnprojectVector(x, full_x);
    full_x += fixed_normal_frames_;   // this is a misnomer, also used for full frames too
    full_x += fixed_interior_frames_;

    // assert
    assert(full_x.size() % 9 == 0);

    return VectorToMatrix(full_x, 9);
}

// Convert the extended frames to variable
Eigen::VectorXd MiNTModel::ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames) {
    // sanity checks
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return Eigen::VectorXd(0);
    }

    if (extended_V_.rows() == 0) {
        std::cerr << "The extended mesh is not set!" << std::endl;
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

// Evaluate the energy
double MiNTModel::ComputeEnergy(const Eigen::MatrixXd& input_frames, Eigen::MatrixXd* virtual_ext_bnd_frames) {
    // Step 0: make sure the mesh has been set
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return 0.0;
    }

    // Step 1: add additional layers for the boundary
    AddAdditionalLayerForBoundary();

    // Step 2: set the selection matrix
    std::vector<bool> free_dofs(extended_mesh_.nTets() * 9, true);
    proj_ = std::make_unique<Projection>(free_dofs);

    // Step 3: The initial guess
    Eigen::MatrixXd ext_frames(extended_mesh_.nTets(), 9);
    ext_frames.setRandom();
    ext_frames.block(0, 0, input_frames.rows(), input_frames.cols()) = input_frames;
    if (virtual_ext_bnd_frames) {
        ext_frames.block(init_mesh_->nTets(), 0, virtual_ext_bnd_frames->rows(), virtual_ext_bnd_frames->cols()) =
            *virtual_ext_bnd_frames;
    } else {
        for (int i = init_mesh_->nTets(); i < extended_mesh_.nTets(); i++) {
            ext_frames.row(i).segment<3>(0) = fixed_normal_frames_.segment<3>(9 * i);
        }
    }

    // Step 4: The objective function
    MiNTEnergy energy_model(extended_V_, extended_mesh_, *init_mesh_);
    energy_model.Precompute();
    energy_model.weights_.w_mint = 1;
    energy_model.weights_.w_unit_norm = 1;
    energy_model.weights_.w_bound = 1;
    energy_model.weights_.w_scaled_jacobian = 1;

    auto obj_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv = nullptr,
                        Eigen::SparseMatrix<double>* hess = nullptr, bool is_proj = false, double global_scale = 1.0,
                        bool recompute = false) {
        Eigen::MatrixXd cur_ext_frames = ConvertVariableToFrames(x);
        std::vector<Eigen::Triplet<double>> hess_T;

        Eigen::VectorXd unit_deriv;
        Eigen::VectorXd integrability_deriv;
        Eigen::VectorXd bound_penalty_deriv;
        Eigen::VectorXd sjac_deriv;

        // energy_model.weights_ = energy_model.weights_.setGlobalScale(global_scale);
        bool is_krush = energy_model.weights_.b_use_kruskal_tensors_for_sym_smoothness;

        double smoothness_energy =
            energy_model.ComputeSmoothnessEnergyWithTriplets(cur_ext_frames, deriv, hess ? &hess_T : nullptr, is_proj);
        double mint_energy = energy_model.ComputeIntegrabilityEnergyWithTriplets(
            cur_ext_frames, deriv ? &integrability_deriv : nullptr, hess ? &hess_T : nullptr, is_proj, is_krush);
        double unit_energy = energy_model.ComputeUnitNormPenaltyWithTriplets(
            cur_ext_frames, deriv ? &unit_deriv : nullptr, hess ? &hess_T : nullptr, is_proj);
        double boundary_penalty_energy = energy_model.ComputeUnitNormalBoundaryPenaltyWithTriplets(
            cur_ext_frames, deriv ? &bound_penalty_deriv : nullptr, hess ? &hess_T : nullptr, is_proj);
        double sjac_energy = energy_model.ComputeScaledJacobianPenaltyWithTriplets(
            cur_ext_frames, deriv ? &sjac_deriv : nullptr, hess ? &hess_T : nullptr, is_proj);

        double energy = 0;
        energy += smoothness_energy;
        energy += mint_energy;
        energy += unit_energy;
        energy += boundary_penalty_energy;
        energy += sjac_energy;

        if (deriv) {
            *deriv += integrability_deriv;
            *deriv += unit_deriv;
            *deriv += bound_penalty_deriv;
            *deriv += sjac_deriv;
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
    Eigen::VectorXd x = ConvertFramesToVariable(ext_frames);

    return obj_func(x, nullptr, nullptr, false, 1.);
}

bool MiNTModel::ShowDensity(SolveParams solve_params, const Eigen::MatrixXd& input_frames,
                            Eigen::MatrixXd& output_frames, Eigen::MatrixXd* virtual_ext_bnd_frames,
                            Eigen::MatrixXd* opt_virtual_ext_bnd_frames) {
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        // Eigen::VectorXd blah;
        return false;
    }

    // Step 1: add additional layers for the boundary
    AddAdditionalLayerForBoundary();

    // Step 2: set the selection matrix
    std::vector<bool> free_dofs(extended_mesh_.nTets() * 9, true);
    for (int i = init_mesh_->nTets(); i < extended_mesh_.nTets(); i++) {
        // free_dofs[9 * i] = false;
        // free_dofs[9 * i + 1] = false;
        // free_dofs[9 * i + 2] = false;
    }
    proj_ = std::make_unique<Projection>(free_dofs);

    // Step 3: The initial guess
    Eigen::MatrixXd ext_frames(extended_mesh_.nTets(), 9);
    ext_frames.setRandom();

    // load the input frames
    Eigen::MatrixXd ext_frames_orig = ext_frames;
    ext_frames.block(0, 0, input_frames.rows(), input_frames.cols()) = input_frames;
    if (virtual_ext_bnd_frames) {
        ext_frames.block(init_mesh_->nTets(), 0, virtual_ext_bnd_frames->rows(), virtual_ext_bnd_frames->cols()) =
            *virtual_ext_bnd_frames;
    } else {
        for (int i = init_mesh_->nTets(); i < extended_mesh_.nTets(); i++) {
            ext_frames.row(i).segment<3>(0) = fixed_normal_frames_.segment<3>(9 * i);
        }
    }

    // Step 4: Define the objective function and corresponding logging lambda
    MiNTEnergy energy_model(extended_V_, extended_mesh_, *init_mesh_);
    energy_model.Precompute();
    energy_model.weights_ = solve_params;
    energy_model.show_energy_ = true;
    energy_model.ComputePermutations(ext_frames);

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeSmoothnessEnergyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_smooth_density = energy_model.per_tet_energy_density_;

    energy_model.weights_.b_smooth_asap_combing = true;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeCombedSmoothnessEnergyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_asap_smoothness_density = energy_model.per_tet_energy_density_;

    // std::cout << cur_asap_smoothness_density.maxCoeff() << std::endl;
    // std::cout << cur_asap_smoothness_density << "asap" << std::endl;

    energy_model.weights_.b_smooth_asap_combing = false;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeCombedSmoothnessEnergyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_aiap_smoothness_density = energy_model.per_tet_energy_density_;
    // std::cout << cur_aiap_smoothness_density << "aiap" << std::endl;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeIntegrabilityEnergyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_mint_density = energy_model.per_tet_energy_density_;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeDeterminantPenaltyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_neohookean_density = energy_model.per_tet_energy_density_;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeUnitNormPenaltyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_unit_norm_density = energy_model.per_tet_energy_density_;

    energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    energy_model.ComputeUnitNormalBoundaryPenaltyWithTriplets(ext_frames, nullptr, nullptr, false);
    cur_bound_align_density = energy_model.per_tet_energy_density_;

    flagged_frames = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    flagged_edges = Eigen::VectorXd::Zero(extended_mesh_.nEdges());

    for (int fid = 0; fid < extended_mesh_.nEdges(); fid++) {
        int tet_id0 = extended_mesh_.faceTet(fid, 0);
        int tet_id1 = extended_mesh_.faceTet(fid, 1);
        if (tet_id0 == -1 || tet_id1 == -1) {
            continue;   // this shouldn't happen I think, but here for robustness.
        }

        // Eigen::VectorXd f = ext_frames.row(extended_mesh_.edgeTet(fid, 0));
        // Eigen::VectorXd g = ext_frames.row(extended_mesh_.edgeTet(fid, 1));

        Eigen::VectorXd f = ext_frames.row(tet_id0);
        Eigen::VectorXd g = ext_frames.row(tet_id1);

        Eigen::MatrixXd face_basis0 = energy_model.tet_facet_basis_.at(tet_id0).at(extended_mesh_.faceTetIndex(fid, 0));
        Eigen::MatrixXd id_basis = Eigen::MatrixXd::Identity(3, 3);

        Eigen::MatrixXd P = ComputeOptimalPermutationMatrix(f, g, face_basis0);
        Eigen::MatrixXd Pi = ComputeOptimalPermutationMatrix(f, g, id_basis);

        Eigen::MatrixXd R = P.transpose() * Pi - Eigen::MatrixXd::Identity(9, 9);

        if (R.norm() > 1e-6) {
            flagged_frames[tet_id0] = 1;
            flagged_frames[tet_id1] = 1;
            flagged_edges[fid] = 1;
        }
    }

    // energy_model.per_tet_energy_density_ = Eigen::VectorXd::Zero(extended_mesh_.nTets());
    // energy_model.ComputeScaledJacobianPenaltyWithTriplets(ext_frames, nullptr, nullptr, false);
    // cur_scaled_jacobian_density = energy_model.per_tet_energy_density_;

    return true;
}

// Solve the smooth frames along with boundary normals
bool MiNTModel::SolveSmoothFrames(SolveParams solve_params, const Eigen::MatrixXd& input_frames,
                                  Eigen::MatrixXd& output_frames, Eigen::MatrixXd* virtual_ext_bnd_frames,
                                  Eigen::MatrixXd* opt_virtual_ext_bnd_frames, std::function<void()>* viz_callback) {
    // Step 0: make sure the mesh has been set
    if (!init_V_ || !init_mesh_) {
        std::cerr << "The initial mesh is not set!" << std::endl;
        return false;
    }

    // Step 1: add additional layers for the boundary
    AddAdditionalLayerForBoundary();

    // override normal boundary conditions with input boundary conditions
    int nInitTets = init_mesh_->nTets();
    fixed_normal_frames_.setZero(9 * (extended_T_.rows()));
    if (virtual_ext_bnd_frames && solve_params.boundary_hard_constraints == BoundaryHardConstraints::Normal) {
        for (int i = 0; i < virtual_ext_bnd_frames->rows(); i++) {
            fixed_normal_frames_.segment<3>(9 * (i + nInitTets)) = virtual_ext_bnd_frames->row(i).segment<3>(0);
        }
    } else if (virtual_ext_bnd_frames && solve_params.boundary_hard_constraints == BoundaryHardConstraints::Frame) {
        for (int i = 0; i < virtual_ext_bnd_frames->rows(); i++) {
            fixed_normal_frames_.segment<9>(9 * (i + nInitTets)) = virtual_ext_bnd_frames->row(i).segment<9>(0);
        }
    }

    // Step 2: set the selection matrix
    std::vector<bool> free_dofs(extended_mesh_.nTets() * 9, true);
    for (int i = init_mesh_->nTets(); i < extended_mesh_.nTets(); i++) {
        if (solve_params.boundary_hard_constraints == BoundaryHardConstraints::Normal) {
            free_dofs[9 * i] = false;
            free_dofs[9 * i + 1] = false;
            free_dofs[9 * i + 2] = false;
        }

        if (solve_params.boundary_hard_constraints == BoundaryHardConstraints::Frame) {
            free_dofs[9 * i] = false;
            free_dofs[9 * i + 1] = false;
            free_dofs[9 * i + 2] = false;
            free_dofs[9 * i + 3] = false;
            free_dofs[9 * i + 4] = false;
            free_dofs[9 * i + 5] = false;
            free_dofs[9 * i + 6] = false;
            free_dofs[9 * i + 7] = false;
            free_dofs[9 * i + 8] = false;
        }
    }

    int rank = 3;
    // setting for constrained interior frames.
    fixed_interior_frames_.setZero(9 * (extended_T_.rows()));
    for (int i = 0; i < input_frames.rows(); i++) {
        // if (rank == 1 || rank == 2) {
        //     fixed_interior_frames_.segment<3>(9 * i + 3) = input_frames.row(i).segment<3>(3);
        // }
        // if (rank == 1)
        // {
        //     fixed_interior_frames_.segment<3>(9 * i + 6) = input_frames.row(i).segment<3>(6);
        // }
    }

    for (int i = 0; i < extended_mesh_.nTets(); i++) {
        if (rank == 1 || rank == 2) {
            free_dofs[9 * i + 3] = false;
            free_dofs[9 * i + 4] = false;
            free_dofs[9 * i + 5] = false;
        }
        if (rank == 1) {
            free_dofs[9 * i + 6] = false;
            free_dofs[9 * i + 7] = false;
            free_dofs[9 * i + 8] = false;
        }
    }
    proj_ = std::make_unique<Projection>(free_dofs);

    // Init Log Folder

    // Store the constructed path in AppState
    logOuterIter = logPath + "/outer_iters";
    logInnerIter = logPath + "/inner_iters";
    std::string dataLogPath = logPath + "/run_data";

    std::string frames_folderpath = "/frames";
    std::string boundary_frames_folderpath = "/boundary_frames";
    std::string state_step_folderpath = "/state_step";

    std::cout << logOuterIter << std::endl;
    std::cout << logInnerIter << std::endl;

    // Create the directory using std::filesystem
    std::filesystem::create_directories(logOuterIter);
    std::filesystem::create_directories(logInnerIter);
    std::filesystem::create_directories(dataLogPath);

    std::filesystem::create_directories(logInnerIter + frames_folderpath);
    std::filesystem::create_directories(logInnerIter + boundary_frames_folderpath);
    std::filesystem::create_directories(logInnerIter + state_step_folderpath);

    // Step 3: Define the objective function and corresponding logging lambda
    MiNTEnergy energy_model(extended_V_, extended_mesh_, *init_mesh_);
    energy_model.Precompute();
    energy_model.weights_ = solve_params;

    // Step 4: The initial guess
    Eigen::MatrixXd ext_frames(extended_mesh_.nTets(), 9);
    ext_frames.setRandom();

    // // Set input boundary frames aligned to crease
    Eigen::MatrixXd input_frames_aligned = input_frames;
    for (int i = 0; i < sharp_feature_edges_.size(); i++) {
        //
        int nconstraints = std::min((int)sharp_feature_edges_[i].size(), 2);
        if (sharp_feature_edges_[i].size() > 2) {
            std::cout << "Warning: More than 2 sharp features detected for boundary tet " << i << std::endl;
            std::cout << "Skipping this constraint" << std::endl;
            std::cout << "sharp_feature_edges_[i] " << sharp_feature_edges_[i].size() << std::endl;
            for (int j = 0; j < sharp_feature_edges_[i].size(); j++) {
                std::cout << sharp_feature_edges_[i][j].transpose() << std::endl;
            }
            std::cout << "consider subdividing the mesh to avoid this issue." << std::endl;
            sharp_feature_edges_[i].resize(0);
        }
        for (int j = 0; j < nconstraints; j++) {
            // std::cout << "i: " << i << "j: " << j << ", nInitTets: " << nInitTets
            //           << ", i + nInitTets: " << i + nInitTets
            //           << ", input_frames_aligned.rows(): " << input_frames_aligned.rows() << std::endl;
            Eigen::Vector3d edge = sharp_feature_edges_[i][j];
            Eigen::Vector3d edge_normal = edge.normalized();
            Eigen::Vector3d b_normals = energy_model.boundary_normals_.row(i);
            Eigen::Vector3d edge_perp = edge_normal.cross(b_normals);

            // std::cout << "edge: " << edge.transpose() << std::endl;
            // std::cout << "edge_normal: " << edge_normal.transpose() << std::endl;
            // std::cout << "b_normals: " << b_normals.transpose() << std::endl;
            // std::cout << "edge_perp: " << edge_perp.transpose() << std::endl;

            assert(3 * j + 3 <= input_frames_aligned.cols());
            assert(i < input_frames_aligned.rows());
            virtual_ext_bnd_frames->row(i).segment<3>(3 * (j + 1)) = edge_perp * 25;
        }
    }

    energy_model.boundary_features_ = sharp_feature_edges_;

    // load the input frames
    Eigen::MatrixXd ext_frames_orig = ext_frames;
    ext_frames.block(0, 0, input_frames.rows(), input_frames.cols()) = input_frames;
    if (virtual_ext_bnd_frames) {
        ext_frames.block(init_mesh_->nTets(), 0, virtual_ext_bnd_frames->rows(), virtual_ext_bnd_frames->cols()) =
            *virtual_ext_bnd_frames;
    } else {
        for (int i = init_mesh_->nTets(); i < extended_mesh_.nTets(); i++) {
            ext_frames.row(i).segment<3>(0) = fixed_normal_frames_.segment<3>(9 * i);
        }
    }

    proj_->setFixedVars(ext_frames);   // apply mask to init frames to obtain hard boundary terms.

    Eigen::MatrixXd reference_field = ext_frames;
    // Eigen::VectorXd idframe = Eigen::VectorXd(9);
    // idframe << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    // for (int i = 0; i < reference_field.rows(); i++) {
    //     reference_field.row(i) = idframe;
    // }
    energy_model.ComputeGuidedTensors(reference_field);

    // This makes it so that boundary gets weighted more on larger meshes
    // so that weights do the same thing under refinement.
    // (double)init_mesh_->nTets() / (double) init_mesh_->nBoundaryElements();
    energy_model.weights_.w_bound_rescale = .075 * std::pow(init_mesh_->nBoundaryElements(), 1.5);

    std::cout << "Solver Input Config: " << std::endl;
    std::cout << "Inner Step Limit: " << energy_model.weights_.inner_iter_max_steps
              << " Grad Norm Exit Tol: " << energy_model.weights_.grad_tol
              << " delta x per step exit tol: " << energy_model.weights_.xTol
              << " delta f per step exit tol: " << energy_model.weights_.fTol << std::endl;

    std::cout << "Energy Weights: " << std::endl;
    std::cout << "mint init: " << energy_model.weights_.w_mint << " mint max: " << energy_model.weights_.w_max_lambda
              << " smooth: " << energy_model.weights_.w_smooth << " unit norm: " << energy_model.weights_.w_unit_norm
              << " unit barrier " << energy_model.weights_.w_unit_barrier
              << " boundary: " << energy_model.weights_.w_bound
              << " scaled jacobian: " << energy_model.weights_.w_scaled_jacobian
              << " boundary rescale: " << energy_model.weights_.w_bound_rescale << std::endl;

    auto run_gradient_hessian_convergence_checks =
        [&](Eigen::MatrixXd& ext_frames, MiNTEnergy& energy_model,
            std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool, double,
                                 bool)>
                obj_func,
            const Eigen::VectorXd& x0) {
            ext_frames = ext_frames + 0.1 * Eigen::MatrixXd::Random(ext_frames.rows(), ext_frames.cols());

            energy_model.frames_prev_outer_step_ =
                ext_frames + 0.01 * Eigen::MatrixXd::Random(ext_frames.rows(), ext_frames.cols());

            energy_model.ComputePermutations(ext_frames);

            // energy_model.TestComputeScaledJacobianPenaltyPerTet(ext_frames, 1);
            energy_model.TestComputeOrthogonalityEnergyPerTet(ext_frames, 1);

            energy_model.TestIntegrabilityMintVsPrimal(ext_frames, true);

            energy_model.TestComputeSmoothnessEnergy(ext_frames);

            energy_model.TestComputeCombedSmoothnessEnergy(ext_frames);
            energy_model.TestComputeCombedIntegrabilityEnergy(ext_frames);

            energy_model.TestComputeIntegrabilityEnergy(ext_frames);

            // energy_model.TestComputeUnitBarrierPenaltyPerTet(ext_frames, 1);

            energy_model.TestComputeUnitNormalBoundaryPenalty(ext_frames);
            energy_model.TestComputeUnitNormPenalty(ext_frames);

            energy_model.TestComputeDeterminantPenaltyPerTet(ext_frames, 1);

            energy_model.TestComputeViscosityPenalty(ext_frames);

            energy_model.TestComputeFeatureAlignmentPenalty(ext_frames);

            // std::cout << x0.transpose() << std::endl;

            OptSolver::TestFuncGradHessian(obj_func, x0);

            return true;
        };

    auto per_inner_iter_update = [&](const Eigen::VectorXd& x) {
        energy_model.frames_prev_inner_step_ = ConvertVariableToFrames(x);
        return 0;
    };

    auto log_obj_func = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd cur_ext_frames = ConvertVariableToFrames(x);

        bool is_proj = false;

        double undo_scale = 1.0 / energy_model.weights_.getGlobalScale();

        double energy = 0;

        // All the primal/combed terms
        // if (!energy_model.weights_.b_headless_mode && false) {
        bool prev_b_smooth_asap_combing = energy_model.weights_.b_smooth_asap_combing;

        // Compute all energy terms
        double primal_integrability_energy = energy_model.ComputePrimalIntegrabilityEnergy(cur_ext_frames);

        energy_model.weights_.b_smooth_asap_combing = true;

        double asap_smoothness_energy =
            energy_model.ComputeCombedSmoothnessEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);

        energy_model.weights_.b_smooth_asap_combing = false;

        double aiap_smoothness_energy =
            energy_model.ComputeCombedSmoothnessEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);

        double scale = energy_model.weights_.getGlobalScale();
        double smoothness_scale =
            (energy_model.weights_.w_smooth_sym * scale + energy_model.weights_.w_smooth_sym_precon * scale * scale);

        asap_smoothness_energy /= smoothness_scale;
        aiap_smoothness_energy /= smoothness_scale;

        energy_model.weights_.b_smooth_asap_combing = prev_b_smooth_asap_combing;

        double combed_dirichlet_energy =
            energy_model.ComputeCombedSmoothnessEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);

        double combed_curl_energy =
            energy_model.ComputeCombedIntegrabilityEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);

        combed_dirichlet_energy *= undo_scale;

        combed_curl_energy *= undo_scale;

        energyLog.get()->asap_combed_smoothness.push_back(asap_smoothness_energy);
        energyLog.get()->aiap_combed_smoothness.push_back(aiap_smoothness_energy);
        energyLog.get()->combed_smoothness.push_back(combed_dirichlet_energy);
        energyLog.get()->combed_integrability.push_back(combed_curl_energy);
        energyLog.get()->primal_integrability.push_back(primal_integrability_energy);

        energyLog.get()->global_scale_outiter.push_back(scale);

        if (energy_model.weights_.b_smooth_combed) energy += combed_dirichlet_energy;
        if (energy_model.weights_.b_int_combed) energy += combed_curl_energy;

        energy_model.weights_.primal_integrability_energy = primal_integrability_energy;
        energy_model.weights_.asap_combed_smoothness_energy = asap_smoothness_energy;
        energy_model.weights_.aiap_combed_smoothness_energy = aiap_smoothness_energy;
        energy_model.weights_.combed_smoothness_energy = combed_dirichlet_energy;
        // }

        double dirichlet_energy =
            energy_model.ComputeSmoothnessEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);
        double mint_energy =
            energy_model.ComputeIntegrabilityEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj, false);
        double unit_energy = energy_model.ComputeUnitNormPenaltyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj);
        double sjac_energy =
            energy_model.ComputeScaledJacobianPenaltyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj);
        double fit_energy = energy_model.ComputeDeviationEnergyWithTriplets(cur_ext_frames, nullptr, nullptr, is_proj);

        // Apply scaling

        dirichlet_energy *= undo_scale;

        mint_energy *= undo_scale;
        unit_energy *= undo_scale;
        sjac_energy *= undo_scale;
        fit_energy *= undo_scale;

        // primal_integrability_energy *= undo_scale;
        // asap_smoothness_energy *= undo_scale;

        // Log energy terms
        energyLog.get()->smoothness.push_back(dirichlet_energy);

        energyLog.get()->symmetric_integrability.push_back(mint_energy);

        energyLog.get()->unit_penalty.push_back(unit_energy);
        energyLog.get()->scaled_jacobian.push_back(sjac_energy);
        energyLog.get()->fit_penalty.push_back(fit_energy);

        // Compute the total energy
        if (energy_model.weights_.b_smooth_sym) energy += dirichlet_energy;
        if (energy_model.weights_.b_int_sym) energy += mint_energy;
        if (energy_model.weights_.b_unit) energy += unit_energy;
        if (energy_model.weights_.b_orthog) energy += sjac_energy;
        if (energy_model.weights_.b_fit) energy += fit_energy;

        // Log the total energy
        energyLog.get()->total_energy_unscale.push_back(energy);

        // Output energy details
        std::cout << "smoothness_energy: " << dirichlet_energy << " mint_energy: " << mint_energy
                  << " unit_energy: " << unit_energy << " orthog_energy: " << sjac_energy << " total energy: " << energy
                  << " undo scale: " << undo_scale << std::endl;

        // Serialize the current frames to file
        currentFileID_++;
        energy_model.weights_.total_step = currentFileID_;
        std::string suffix = std::to_string(currentFileID_ + 100000);

        Eigen::MatrixXd log_frames = cur_ext_frames.block(0, 0, input_frames.rows(), 9);
        Eigen::MatrixXd log_boundary_frames =
            cur_ext_frames.block(input_frames.rows(), 0, (extended_mesh_.nTets() - init_mesh_->nTets()), 9);

        Serialization::serializeMatrix(log_frames, logInnerIter + "/frames/f" + "_" + suffix + ".bfra", 3);
        Serialization::serializeMatrix(log_boundary_frames,
                                       logInnerIter + "/boundary_frames/bf" + "_" + suffix + ".bfra", 3);

        // Update and save the energy model weights
        energy_model.weights_.smoothness_energy = dirichlet_energy;
        energy_model.weights_.mint_energy = mint_energy;
        energy_model.weights_.unit_energy = unit_energy;

        energy_model.weights_.scaled_jacobian_energy = sjac_energy;
        energy_model.weights_.total_energy = energy;

        Serialization::serializeSolveParams(energy_model.weights_,
                                            logInnerIter + "/state_step/conf" + "_" + suffix + ".json");

        return energy;
    };

    auto obj_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv = nullptr,
                        Eigen::SparseMatrix<double>* hess = nullptr, bool is_proj = false, double global_scale = 1.,
                        bool recompute_state = false) {
        Eigen::MatrixXd cur_ext_frames = ConvertVariableToFrames(x);
        std::vector<Eigen::Triplet<double>> hess_T;

        std::vector<Eigen::Triplet<double>> unit_hess;
        std::vector<Eigen::Triplet<double>> mint_hess;

        std::vector<Eigen::Triplet<double>> dirichlet_hess;
        std::vector<Eigen::Triplet<double>> bound_hess;
        std::vector<Eigen::Triplet<double>> sjac_hess;
        std::vector<Eigen::Triplet<double>> fit_hess;
        std::vector<Eigen::Triplet<double>> visc_hess;
        std::vector<Eigen::Triplet<double>> feat_align_hess;

        std::vector<Eigen::Triplet<double>> combed_dirichlet_hess;
        std::vector<Eigen::Triplet<double>> combed_div_hess;
        std::vector<Eigen::Triplet<double>> combed_curl_hess;

        Eigen::VectorXd unit_deriv;
        Eigen::VectorXd mint_deriv;
        Eigen::VectorXd dirichlet_deriv;

        Eigen::VectorXd combed_dirichlet_deriv;
        Eigen::VectorXd combed_div_deriv;
        Eigen::VectorXd combed_curl_deriv;

        Eigen::VectorXd sjac_deriv;
        Eigen::VectorXd det_deriv;
        Eigen::VectorXd fit_deriv;
        Eigen::VectorXd visc_deriv;
        Eigen::VectorXd feat_align_deriv;

        if (recompute_state) {
            energy_model.ComputePermutations(cur_ext_frames);
        }

        bool is_krush = energy_model.weights_.b_use_kruskal_tensors_for_sym_smoothness;

        // energy_model.weights_ = energy_model.weights_.setGlobalScale(global_scale);

        double mint_energy;
        double dirichlet_energy;
        double combed_dirichlet_energy;
        double combed_curl_energy;
        double unit_energy;
        double sjac_energy;
        // double det_energy;
        double fit_energy;
        double visc_energy;
        double feat_align_energy;

        if (energy_model.weights_.b_smooth_sym) {
            dirichlet_energy = energy_model.ComputeSmoothnessEnergyWithTriplets(
                cur_ext_frames, deriv ? &dirichlet_deriv : nullptr, hess ? &dirichlet_hess : nullptr, is_proj, false);
        }

        if (energy_model.weights_.b_smooth_combed) {
            combed_dirichlet_energy = energy_model.ComputeCombedSmoothnessEnergyWithTriplets(
                cur_ext_frames, deriv ? &combed_dirichlet_deriv : nullptr, hess ? &combed_dirichlet_hess : nullptr,
                is_proj, false);
        }

        if (energy_model.weights_.b_int_sym) {
            mint_energy = energy_model.ComputeIntegrabilityEnergyWithTriplets(
                cur_ext_frames, deriv ? &mint_deriv : nullptr, hess ? &mint_hess : nullptr, is_proj, false);
        }

        if (energy_model.weights_.b_int_combed) {
            combed_curl_energy = energy_model.ComputeCombedIntegrabilityEnergyWithTriplets(
                cur_ext_frames, deriv ? &combed_curl_deriv : nullptr, hess ? &combed_curl_hess : nullptr, is_proj,
                false);
        }

        if (energy_model.weights_.b_unit) {
            unit_energy = energy_model.ComputeUnitNormPenaltyWithTriplets(cur_ext_frames, deriv ? &unit_deriv : nullptr,
                                                                          hess ? &unit_hess : nullptr, is_proj);
        }

        if (energy_model.weights_.b_orthog) {
            sjac_energy = energy_model.ComputeScaledJacobianPenaltyWithTriplets(
                cur_ext_frames, deriv ? &sjac_deriv : nullptr, hess ? &sjac_hess : nullptr, is_proj);
        }

        if (energy_model.weights_.b_fit) {
            fit_energy = energy_model.ComputeDeviationEnergyWithTriplets(cur_ext_frames, deriv ? &fit_deriv : nullptr,
                                                                         hess ? &fit_hess : nullptr, is_proj);
            // det_energy = energy_model.ComputeDeterminantPenaltyWithTriplets(cur_ext_frames, deriv ? &det_deriv :
            // nullptr, hess ? &det_hess : nullptr, is_proj);
        }

        if (energy_model.weights_.b_viscosity) {
            visc_energy = energy_model.ComputeViscosityPenaltyWithTriplets(
                cur_ext_frames, deriv ? &visc_deriv : nullptr, hess ? &visc_hess : nullptr, is_proj);
        }

        if (energy_model.weights_.b_feat_align) {
            feat_align_energy = energy_model.ComputeFeatureAlignmentPenaltyWithTriplets(
                cur_ext_frames, deriv ? &feat_align_deriv : nullptr, hess ? &feat_align_hess : nullptr, is_proj);
        }

        double energy = 0;

        if (energy_model.weights_.b_smooth_sym) {
            energy += dirichlet_energy;
        }
        if (energy_model.weights_.b_smooth_combed) {
            energy += combed_dirichlet_energy;
        }
        if (energy_model.weights_.b_int_sym) {
            energy += mint_energy;
        }
        if (energy_model.weights_.b_int_combed) {
            energy += combed_curl_energy;
        }
        if (energy_model.weights_.b_unit) {
            energy += unit_energy;
        }
        if (energy_model.weights_.b_orthog) {
            energy += sjac_energy;
        }
        if (energy_model.weights_.b_fit) {
            energy += fit_energy;
        }
        if (energy_model.weights_.b_viscosity) {
            energy += visc_energy;
        }
        if (energy_model.weights_.b_feat_align) {
            energy += feat_align_energy;
        }

        if (deriv) {
            deriv->setZero(cur_ext_frames.rows() * cur_ext_frames.cols());

            if (energy_model.weights_.b_smooth_sym) {
                *deriv += dirichlet_deriv;
            }
            if (energy_model.weights_.b_smooth_combed) {
                *deriv += combed_dirichlet_deriv;
            }
            if (energy_model.weights_.b_int_sym) {
                *deriv += mint_deriv;
            }
            if (energy_model.weights_.b_int_combed) {
                *deriv += combed_curl_deriv;
            }
            if (energy_model.weights_.b_unit) {
                *deriv += unit_deriv;
            }
            if (energy_model.weights_.b_orthog) {
                *deriv += sjac_deriv;
            }
            if (energy_model.weights_.b_fit) {
                *deriv += fit_deriv;
            }
            if (energy_model.weights_.b_viscosity) {
                *deriv += visc_deriv;
            }
            if (energy_model.weights_.b_feat_align) {
                *deriv += feat_align_deriv;
            }

            Eigen::VectorXd proj_deriv;
            proj_->ProjectVector(*deriv, proj_deriv);
            *deriv = std::move(proj_deriv);
        }

        if (hess) {
            double hsize = mint_hess.size() + combed_dirichlet_hess.size() + dirichlet_hess.size() +
                           combed_curl_hess.size() + unit_hess.size() + sjac_hess.size();
            hess_T.reserve(hsize);
            hess_T.insert(hess_T.end(), unit_hess.begin(), unit_hess.end());
            hess_T.insert(hess_T.end(), mint_hess.begin(), mint_hess.end());
            hess_T.insert(hess_T.end(), dirichlet_hess.begin(), dirichlet_hess.end());
            hess_T.insert(hess_T.end(), combed_dirichlet_hess.begin(), combed_dirichlet_hess.end());
            hess_T.insert(hess_T.end(), combed_curl_hess.begin(), combed_curl_hess.end());
            hess_T.insert(hess_T.end(), sjac_hess.begin(), sjac_hess.end());
            hess_T.insert(hess_T.end(), fit_hess.begin(), fit_hess.end());
            hess_T.insert(hess_T.end(), visc_hess.begin(), visc_hess.end());
            hess_T.insert(hess_T.end(), feat_align_hess.begin(), feat_align_hess.end());
            proj_->ProjectMatrix(hess_T);
            int nvars = proj_->ProjDOFs();
            hess->resize(nvars, nvars);
            hess->setFromTriplets(hess_T.begin(), hess_T.end());
        }

        return energy;
    };

    //  auto obj_func = normalAlignOpt(energy_model);

    // Step 5: Find the max step size function
    auto find_max_step = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& dir) {
        return 1.0;
    };

    // Step 5: The optimization
    Eigen::VectorXd x0 = ConvertFramesToVariable(ext_frames);
    energy_model.frames_prev_outer_step_ = ext_frames;
    energy_model.frames_prev_inner_step_ = ext_frames;
    bool to_proj = true;

    energy_model.ComputePermutations(ext_frames);

    if (energy_model.weights_.w_fit * energy_model.weights_.w_self_align > 1e-10) {
        std::cout << "error, not supported fitting and using self align energy because of a hack" << std::endl;
        std::cout << "if w_fit > 0, fit to init frame, if w_self_align > 0 fit to each outer iter" << std::endl;

        return false;
    }

    // Test_ComputeOptimalPermutationMatrix();

    // return false;

    // run_gradient_hessian_convergence_checks(ext_frames, energy_model, obj_func, ConvertFramesToVariable(ext_frames));

    // return false;

    // Actual optimization

    int outer_iter = 0;

    double w_mint_init = energy_model.weights_.w_int_sym;

    // so that penalty stops at max_lambda
    // this actual lambda paramter might exceed the max which is confusing.
    energy_model.weights_.w_max_lambda *= 1. / w_mint_init;

    double min_scale = 1. / energy_model.weights_.w_max_lambda;

    // energy_model.weights_.w_mint = 1e-8;

    double scale = 1;
    double no_int_converged = 0;

    bool is_initalizing = true;
    // int max_interior_iters = energy_model.weights_.inner_iter_max_steps;
    // energy_model.weights_.inner_iter_max_steps = 50;

    ///////
    //
    // SOLVER LOOP
    //
    ///////
    log_obj_func(x0);
    while (energy_model.weights_.lambda_penalty <
           energy_model.weights_.w_max_lambda)   // 1e10)   // 1e0) //  1e8)  one_step 1e6 1e20
    {
        // if ( energy_model.weights_.w_mint < .999 )
        // {
        double mult_factor = energy_model.weights_.w_outer_step;
        energy_model.weights_.lambda_penalty *= mult_factor;

        if (energy_model.weights_.lambda_penalty * w_mint_init < .01) {
            // energy_model.weights_.lambda_penalty *= 10.;
        } else {
            is_initalizing = false;
            // energy_model.weights_.inner_iter_max_steps = max_interior_iters;
        }

        energy_model.weights_.lambda_penalty *= mult_factor;
        double lambda_penalty = energy_model.weights_.lambda_penalty;
        energy_model.weights_.w_mint = w_mint_init * lambda_penalty * energy_model.weights_.getGlobalScale();   // 100

        energy_model.weights_ = energy_model.weights_.setGlobalScale(1. / (w_mint_init * lambda_penalty));

        // this sets the visc term to the solution of the prev outer step
        // energy_model.frames_prev_inner_step_ =  energy_model.frames_prev_outer_step_;
        energy_model.frames_prev_outer_step_ = ConvertVariableToFrames(x0);

        // HACK FLAG THIS OFF!!!!!
        if (energy_model.weights_.b_self_align) {
            energy_model.ComputeGuidedTensors(energy_model.frames_prev_outer_step_);

            // energy_model.weights_.b_self_align = false;
            // std::cout << "error, not supported fitting and using self align energy because of a hack" << std::endl;
            // return false;
        }

        // always start outer iter with all frames with positive volume.
        bool has_inverted_frames = false;
        for (int i = 0; i < energy_model.frames_prev_outer_step_.rows(); i++) {
            Eigen::Vector3d f0 = energy_model.frames_prev_outer_step_.row(i).segment<3>(0);
            Eigen::Vector3d f1 = energy_model.frames_prev_outer_step_.row(i).segment<3>(3);
            Eigen::Vector3d f2 = energy_model.frames_prev_outer_step_.row(i).segment<3>(6);

            double det = (f0.cross(f1)).dot(f2);

            if (det < 0) {
                energy_model.frames_prev_outer_step_.row(i).segment<3>(3) = -f1;
                has_inverted_frames = true;
                // x0.segment<3>(9 * i + 3) = -x0.segment<3>(9 * i + 3); // not compatible with the hard constraint
                // formulation
            }
        }
        if (has_inverted_frames) {
            std::cout << "Inverted frames detected.  Correcting." << std::endl;
            x0 = ConvertFramesToVariable(energy_model.frames_prev_outer_step_);
        }

        std::cout << "**************************" << std::endl
                  << " OUTER STEP: " << energy_model.weights_.total_step
                  << " Iterate w_mint. NOW SET TO: " << energy_model.weights_.w_mint
                  << " increment set to: " << energy_model.weights_.w_outer_step
                  << " current global scale: " << energy_model.weights_.getGlobalScale() << std::endl
                  << "**************************" << std::endl;
        // log_obj_func(x0);

        // a bit of a hack... can do something more principled eventually...
        energy_model.weights_.inner_iter_max_steps = 25;
        if (energy_model.weights_.getGlobalScale() < 100. && energy_model.weights_.getGlobalScale() > 1e-10) {
            energy_model.weights_.inner_iter_max_steps = 100;
        }

        bool is_success = OptSolver::NewtonSolver(obj_func, log_obj_func, find_max_step, per_inner_iter_update, x0,
                                                  energy_model.weights_.reg, energy_model.weights_.inner_iter_max_steps,
                                                  energy_model.weights_.getGlobalScale(),
                                                  energy_model.weights_.grad_tol, energy_model.weights_.xTol,
                                                  energy_model.weights_.fTol, to_proj, true, false, solveData.get());

        solveData.get()->logToFile(dataLogPath);
        energyLog.get()->logToFile(dataLogPath);

        if (!is_success) {
            std::cout << "Newton Solver Failed. Exiting." << std::endl
                      << " current integrability weight: " << energy_model.weights_.w_mint << std::endl;
            break;
        }
        // to_proj = true;
        energy_model.weights_.outer_step++;

        (*viz_callback)();

        // log_obj_func(x0);

        if (solveData.get()->energy_diff.back() < 1e-20 && energy_model.weights_.w_mint < 1e-20) {
            no_int_converged++;
            if (no_int_converged > 10) {
                std::cout << "vanishing energy diff, breaking out " << std::endl;
                break;
            }
        }

        outer_iter++;
        // log less in headless mode if necessary.
        // if (outer_iter % 8 == 0 || is_initalizing || !solve_params.b_headless_mode) {
        log_obj_func(x0);
        // }

        // std::cout << energy_model.frames_prev_outer_step_ << std::endl;
    }

    log_obj_func(x0);

    // Step 6: Return
    ext_frames = ConvertVariableToFrames(x0);
    // ext_frames = ext_frames_orig;
    output_frames = ext_frames.block(0, 0, input_frames.rows(), input_frames.cols());

    if (opt_virtual_ext_bnd_frames) {
        *opt_virtual_ext_bnd_frames = ext_frames.block(
            init_mesh_->nTets(), 0, extended_mesh_.nTets() - init_mesh_->nTets(), output_frames.cols());
    }

    std::cout << "blah blah " << std::endl;
    // std::cout << *opt_virtual_ext_bnd_frames << std::endl;

    return true;
}

}   // namespace MiNT3D
