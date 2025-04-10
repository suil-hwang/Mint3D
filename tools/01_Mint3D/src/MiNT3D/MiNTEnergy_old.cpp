#include "MiNTEnergy.h"

#include <fstream>
#include <iostream>

#include <tbb/parallel_for.h>

#include "MiNTCommonFunc.h"

#include "../SymmetricTensor/SymmetricKruskalTensor.h"
#include "../SymmetricTensor/SymmetricMomentTensor.h"

#include <unsupported/Eigen/KroneckerProduct>

namespace MiNT3D {
/*

// precompute the data structures for the energy computation, including the dual volumes of the faces, and dual vertex
// positions (currently we use the volume centroid as the dual vertex position)
void MiNTEnergy::Precompute() {
    // compute the volumes of the tets, and dual vertex positions
    tet_volumes_.resize(mesh_->nTets(), 0.0);
    Eigen::MatrixXd dual_vertex_pos_(mesh_->nTets(), 3);

    for (int i = 0; i < mesh_->nTets(); i++) {
        Eigen::Vector3d v0 = V_->row(mesh_->tetVertex(i, 0));
        Eigen::Vector3d v1 = V_->row(mesh_->tetVertex(i, 1));
        Eigen::Vector3d v2 = V_->row(mesh_->tetVertex(i, 2));
        Eigen::Vector3d v3 = V_->row(mesh_->tetVertex(i, 3));
        tet_volumes_[i] = std::abs((v1 - v0).cross(v2 - v0).dot(v3 - v0)) / 6.0;

        // we now just use the volume centroid as the dual vertex position
        // the circumcenter may be a better choice
        dual_vertex_pos_.row(i) = (v0 + v1 + v2 + v3).transpose() / 4.0;
    }

    // compute the dual volumes of the faces
    dual_face_volumes_.resize(mesh_->nFaces(), 0);
    dual_edge_lengths_.resize(mesh_->nFaces(), 0);
    for (int i = 0; i < mesh_->nFaces(); i++) {
        bool is_bnd_face = false;
        for (int j = 0; j < 2; j++) {
            int tet_id = mesh_->faceTet(i, j);
            if (tet_id == -1) {
                is_bnd_face = true;
                continue;
            } else {
                dual_face_volumes_[i] += tet_volumes_[tet_id] / 4.0;
            }
        }
        if (!is_bnd_face) {
            dual_face_volumes_[i] /= 2.0;

            // compute the dual edge lengths
            dual_edge_lengths_[i] =
                (dual_vertex_pos_.row(mesh_->faceVertex(i, 0)) - dual_vertex_pos_.row(mesh_->faceVertex(i, 1))).norm();
        }
    }

    ave_edge_len_ = 0;
    for (int i = 0; i < mesh_->nEdges(); i++) {
        int vid0 = mesh_->edgeVertex(i, 0);
        int vid1 = mesh_->edgeVertex(i, 1);
        ave_edge_len_ += (V_->row(vid0) - V_->row(vid1)).norm();
    }
    ave_edge_len_ /= mesh_->nEdges();

    // we need 2nd-, 4th-, 6th-order tensors for 3 dimensional vectors
    tensor_ord_weight_[2] = SymmetricTensor::Get2ndOrderCoeffWeights(3);
    tensor_ord_weight_[4] = SymmetricTensor::Get4thOrderCoeffWeights(3);
    tensor_ord_weight_[6] = SymmetricTensor::Get6thOrderCoeffWeights(3);

    // we need 2nd-, 4th-, 6th-order tensors for 3 dimensional vectors
    tensor2D_ord_weight_[2] = SymmetricTensor::Get2ndOrderCoeffWeights(2);
    tensor2D_ord_weight_[4] = SymmetricTensor::Get4thOrderCoeffWeights(2);
    tensor2D_ord_weight_[6] = SymmetricTensor::Get6thOrderCoeffWeights(2);

    // Compute orthogonal per tet facet basis

    ComputePerTetFacetBasis();
    ComputeBoundaryNormals();
}

void MiNTEnergy::ComputePerTetFacetBasis() {
    int ntets = mesh_->nTets();
    tet_facet_basis_.resize(ntets);
    tet_facet_dual_basis_.resize(ntets);
    tet_facet_full_basis_.resize(ntets);

    for (int tetidx = 0; tetidx < ntets; tetidx++) {
        std::vector<Eigen::MatrixXd> cur_tet_bases(4);
        std::vector<Eigen::MatrixXd> cur_tet_dual_bases(4);
        std::vector<Eigen::MatrixXd> cur_tet_full_bases(4);

        for (int idx = 0; idx < 4; idx++) {
            int face_idx = mesh_->tetFace(tetidx, idx);

            Eigen::Matrix3d rot_facet_to_template;
            Eigen::MatrixXd cur_facet_full_basis = Eigen::MatrixXd::Zero(3, 3);

            Eigen::Vector3d a = V_->row(mesh_->faceVertex(face_idx, 0));
            Eigen::Vector3d b = V_->row(mesh_->faceVertex(face_idx, 1));
            Eigen::Vector3d c = V_->row(mesh_->faceVertex(face_idx, 2));

            Eigen::Vector3d b1 = (b - a).normalized();
            Eigen::Vector3d b2 = (c - a).normalized();
            Eigen::Vector3d n = b1.cross(b2).normalized();

            cur_facet_full_basis.row(0) = b1;
            cur_facet_full_basis.row(1) = b1.cross(n);   // b2
            cur_facet_full_basis.row(2) = n;

            cur_tet_bases[idx] = cur_facet_full_basis.block<2, 3>(0, 0);
            cur_tet_dual_bases[idx] = n.transpose();
            cur_tet_full_bases[idx] = cur_facet_full_basis;
        }
        tet_facet_basis_[tetidx] = cur_tet_bases;
        tet_facet_dual_basis_[tetidx] = cur_tet_dual_bases;
        tet_facet_full_basis_[tetidx] = cur_tet_full_bases;
    }
}

void MiNTEnergy::ComputeBoundaryNormals() {
    boundary_normals_ = Eigen::MatrixXd::Zero(mesh_not_extended_->nBoundaryElements(), 3);
    boundary_b1_ = Eigen::MatrixXd::Zero(mesh_not_extended_->nBoundaryElements(), 3);
    boundary_b2_ = Eigen::MatrixXd::Zero(mesh_not_extended_->nBoundaryElements(), 3);

    for (int i = 0; i < boundary_normals_.rows(); i++) {
        int boundaryFace = mesh_not_extended_->boundaryFace(i);
        Eigen::Vector3d a = V_->row(mesh_not_extended_->faceVertex(boundaryFace, 0));
        Eigen::Vector3d b = V_->row(mesh_not_extended_->faceVertex(boundaryFace, 1));
        Eigen::Vector3d c = V_->row(mesh_not_extended_->faceVertex(boundaryFace, 2));

        Eigen::Vector3d b1 = (b - a).normalized();
        Eigen::Vector3d b2 = (c - a).normalized();
        Eigen::Vector3d n = b1.cross(b2).normalized();

        // appState->bound_centroids.row(i) = ( a + b + c )  / 3.0;
        boundary_normals_.row(i) = n;
        boundary_b1_.row(i) = b1;
        boundary_b2_.row(i) = n.cross(b1);
    }
}

// ComputePermutations(const Eigen::MatrixXd& frames);
void MiNTEnergy::ComputePermutations(const Eigen::MatrixXd& frames) {
    assert(mesh_ && V_);

    permutation_as_smooth_as_possible_.resize(mesh_->nFaces());
    permutation_as_integrable_as_possible_.resize(mesh_->nFaces());

    Eigen::MatrixXd id_basis = Eigen::MatrixXd::Identity(3, 3);

    int nfaces = mesh_->nFaces();

    tbb::parallel_for(tbb::blocked_range<int>(0u, nfaces), [&](const tbb::blocked_range<int>& range) {
        for (int face_id = range.begin(); face_id != range.end(); ++face_id) {
            int tet_id0 = mesh_->faceTet(face_id, 0);
            int tet_id1 = mesh_->faceTet(face_id, 1);
            if (tet_id0 == -1 || tet_id1 == -1) {
                continue;   // Boundary tetrahedron
            }

            Eigen::VectorXd v = frames.row(tet_id0);
            Eigen::VectorXd u = frames.row(tet_id1);

            Eigen::MatrixXd face_basis0 = tet_facet_basis_.at(tet_id0).at(mesh_->faceTetIndex(face_id, 0));

            permutation_as_smooth_as_possible_[face_id] =
                ComputeOptimalPermutationMatrix(frames.row(tet_id0), frames.row(tet_id1), id_basis);

            permutation_as_integrable_as_possible_[face_id] =
                ComputeOptimalPermutationMatrix(frames.row(tet_id0), frames.row(tet_id1), face_basis0);
        }
    });
}

// Compute Guided Tensors
void MiNTEnergy::ComputeGuidedTensors(const Eigen::MatrixXd& guided_frames) {
    assert(mesh_ && V_);
    for (int k = 2; k <= 6; k += 2) {
        // guided_tensor_coeffs_[k] = ComputeNormalizedMomentTensors(guided_frames, k);

        switch (this->weights_.fit_type) {
            case FitType::Moments:

                guided_tensor_coeffs_[k] = ComputeMomentTensors(guided_frames, k);
                // guided_tensor_coeffs_[k] = ComputedMomentTensors(guided_frames, k);
                // throw std::runtime_error("Direction fitting not implemented here yet");

                break;
            case FitType::Direction:
                guided_tensor_coeffs_[k] = ComputeNormalizedMomentTensors(guided_frames, k);
                // tensor = ComputeNormalizedMomentTensorsPerTet(
                //     frames, tet_id, order, (deriv || hess) ? &tensor_deriv : nullptr, hess ? &tensor_hess : nullptr);
                break;
            case FitType::Kruskal:
                guided_tensor_coeffs_[k] = ComputeKruskalTensors(guided_frames, k);
                // ComputeKruskalTensors(const Eigen::MatrixXd& frames, int order,
                //                                                    std::vector<Eigen::MatrixXd>* deriv,
                //                                                    std::vector<std::vector<Eigen::MatrixXd>>* hess)

                // SymmetricTensor::GetMomentCoefficients(frames.row(tet_id), order, deriv ? &tensor_deriv : nullptr,
                //                                                           hess ? &tensor_hess : nullptr);
                // throw std::runtime_error("Kruskal tensor is not implemented here yet");
                break;
            default:
                break;
        }

        for (const auto& matrix : guided_tensor_coeffs_[k]) {
            // std::cout << matrix << std::endl;
        }
        // std::cout << guided_tensor_coeffs_[k] << std::endl;
    }
}

// Compute the Moment tensors
std::vector<Eigen::VectorXd> MiNTEnergy::ComputeMomentTensors(const Eigen::MatrixXd& frames, int order,
                                                              std::vector<Eigen::MatrixXd>* deriv,
                                                              std::vector<std::vector<Eigen::MatrixXd>>* hess) {
    int ntets = mesh_->nTets();

    std::vector<Eigen::VectorXd> tensors(ntets);
    if (deriv) {
        deriv->resize(ntets);
    }
    if (hess) {
        hess->resize(ntets);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            tensors[i] = ComputeMomentTensorsPerTet(frames, i, order, deriv ? &(*deriv)[i] : nullptr,
                                                    hess ? &(*hess)[i] : nullptr);
        }
    });

    return tensors;
}

// Compute the Kruskal tensor of a frame field on a tetrahedron
Eigen::VectorXd MiNTEnergy::ComputeMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                       Eigen::MatrixXd* deriv, std::vector<Eigen::MatrixXd>* hess) {
    int dim = 3;
    assert(frames.cols() % 3 == 0);
    int frame_per_tet = frames.cols() / dim;
    assert(frame_per_tet == 3);

    int coeff_size = SymmetricTensor::GetCoeffSize(dim, order);

    Eigen::VectorXd moment_tensor = Eigen::VectorXd::Zero(coeff_size);
    if (deriv) {
        deriv->resize(coeff_size, frames.cols());
        deriv->setZero();
    }

    if (hess) {
        hess->resize(coeff_size, Eigen::MatrixXd::Zero(frames.cols(), frames.cols()));
    }

    for (int i = 0; i < frame_per_tet; i++) {
        Eigen::VectorXd vec = frames.row(tet_id).segment(dim * i, dim).transpose();
        if (vec.norm() < 1e-12) {
            continue;
        }

        Eigen::MatrixXd deriv_i;
        std::vector<Eigen::MatrixXd> hess_i;
        Eigen::VectorXd moment_tensor_i =
            SymmetricTensor::GetMomentCoefficients(vec, order, deriv ? &deriv_i : nullptr, hess ? &hess_i : nullptr);

        moment_tensor += moment_tensor_i;
        if (deriv) {
            deriv->block(0, dim * i, coeff_size, dim) += deriv_i;
        }
        if (hess) {
            for (int j = 0; j < hess_i.size(); j++) {
                (*hess)[j].block(dim * i, dim * i, dim, dim) += hess_i[j];
            }
        }
    }

    return moment_tensor;
}

// Compute the Kruskal tensors
std::vector<Eigen::VectorXd> MiNTEnergy::ComputeKruskalTensors(const Eigen::MatrixXd& frames, int order,
                                                               std::vector<Eigen::MatrixXd>* deriv,
                                                               std::vector<std::vector<Eigen::MatrixXd>>* hess) {
    int ntets = mesh_->nTets();

    std::vector<Eigen::VectorXd> tensors(ntets);
    if (deriv) {
        deriv->resize(ntets);
    }
    if (hess) {
        hess->resize(ntets);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            tensors[i] = ComputeKruskalTensorsPerTet(frames, i, order, deriv ? &(*deriv)[i] : nullptr,
                                                     hess ? &(*hess)[i] : nullptr);
        }
    });

    // this can be parallelized
    // for (int i = 0; i < ntets; i++) {
    //         tensors[i] = ComputeKruskalTensorsPerTet(frames, i, order, deriv ? &(*deriv)[i] : nullptr, hess ?
    // &(*hess)[i] : nullptr);
    // }

    return tensors;
}

// Compute the Kruskal tensor of a frame field on a tetrahedron
Eigen::VectorXd MiNTEnergy::ComputeKruskalTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                        Eigen::MatrixXd* deriv, std::vector<Eigen::MatrixXd>* hess) {
    int dim = 3;
    assert(frames.cols() % 3 == 0);
    int frame_per_tet = frames.cols() / dim;
    assert(frame_per_tet == 3);

    int coeff_size = SymmetricTensor::GetCoeffSize(dim, order);

    Eigen::VectorXd kruskal_tensor = Eigen::VectorXd::Zero(coeff_size);
    if (deriv) {
        deriv->resize(coeff_size, frames.cols());
        deriv->setZero();
    }

    if (hess) {
        hess->resize(coeff_size, Eigen::MatrixXd::Zero(frames.cols(), frames.cols()));
    }

    for (int i = 0; i < frame_per_tet; i++) {
        Eigen::VectorXd vec = frames.row(tet_id).segment(dim * i, dim).transpose();
        if (vec.norm() < 1e-12) {
            continue;
        }

        Eigen::MatrixXd deriv_i;
        std::vector<Eigen::MatrixXd> hess_i;
        Eigen::VectorXd kruskal_tensor_i =
            SymmetricTensor::GetKruskalCoefficients(vec, order, deriv ? &deriv_i : nullptr, hess ? &hess_i : nullptr);

        kruskal_tensor += kruskal_tensor_i;
        if (deriv) {
            deriv->block(0, dim * i, coeff_size, dim) += deriv_i;
        }
        if (hess) {
            for (int j = 0; j < hess_i.size(); j++) {
                (*hess)[j].block(dim * i, dim * i, dim, dim) += hess_i[j];
            }
        }
    }

    return kruskal_tensor;
}

// Compute the frames projected to each tet face and then lifted as krushkal tensors
std::vector<Eigen::VectorXd> MiNTEnergy::ComputeFaceProjectedTensors(const Eigen::MatrixXd& frames, bool is_kruskal,
                                                                     int order, std::vector<Eigen::MatrixXd>* deriv,
                                                                     std::vector<std::vector<Eigen::MatrixXd>>* hess) {
    int ntets = mesh_->nTets();

    std::vector<Eigen::VectorXd> tensors(ntets * 4);
    if (deriv) {
        deriv->resize(ntets * 4);
    }
    if (hess) {
        hess->resize(ntets * 4);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            for (int j = 0; j < 4; j++) {
                tensors[i * 4 + j] = ComputeBasisProjectedTensorsPerTetFace(
                    tet_facet_basis_.at(i).at(j), frames, is_kruskal, i, order, deriv ? &(*deriv)[i * 4 + j] : nullptr,
                    hess ? &(*hess)[i * 4 + j] : nullptr);
            }
        }
    });

    // for (int i = 0; i != ntets; ++i) {
    // 	for(int j = 0; j < 4; j++)
    // 	{
    // 		tensors[i*4 + j] = ComputeBasisProjectedTensorsPerTetFace(tet_facet_basis_.at(i).at(j), frames,
    // is_kruskal, i, order, deriv ? &(*deriv)[i*4 + j] : nullptr, hess ? &(*hess)[i*4 + j] : nullptr);
    // 	}
    // }

    // if (hess) {
    //     std::cout << "tensors" << std::endl;
    //     for (int fid = 0; fid < mesh_->nFaces(); fid++) {
    //         int tet_id0 = mesh_->faceTet(fid, 0);
    //         int tet_id1 = mesh_->faceTet(fid, 1);
    //         if (tet_id0 == -1 || tet_id1 == -1) {
    //             continue;
    //         }
    //         int idx0 = mesh_->faceTet(fid, 0) * 4 + mesh_->faceTetIndex(fid,0);
    //         int idx1 = mesh_->faceTet(fid, 1) * 4 + mesh_->faceTetIndex(fid,1);
    //         std::cout << "idx0: " << idx0 << " t: " << tensors.at(idx0).transpose() << " | idx1: " << idx1 << " t: "
    //         << tensors.at(idx1).transpose() << std::endl;

    //         std::cout << tet_facet_basis_.at( mesh_->faceTet(fid, 0) ).at(mesh_->faceTetIndex(fid,0)) << std::endl;
    //         std::cout << tet_facet_basis_.at( mesh_->faceTet(fid, 1) ).at(mesh_->faceTetIndex(fid,1)) << std::endl;
    //         std::cout << frames.row(mesh_->faceTet(fid, 0)) << std::endl;
    //         std::cout << frames.row(mesh_->faceTet(fid, 1)) << std::endl;

    //     }
    //     std::cout << "tensors" << std::endl;
    // }

    return tensors;
}

// // Compute the frames projected to each tet face and then lifted as krushkal tensors
// std::vector<Eigen::VectorXd> MiNTEnergy::ComputeFaceProjectedTensors(const Eigen::MatrixXd& frames, bool is_kruskal,
//                                                                      int order, std::vector<Eigen::MatrixXd>* deriv,
//                                                                      std::vector<std::vector<Eigen::MatrixXd>>* hess)
//                                                                      {
//     int ntets = mesh_->nTets();

//     std::vector<Eigen::VectorXd> tensors(ntets * 4);
//     if (deriv) {
//         deriv->resize(ntets * 4);
//     }
//     if (hess) {
//         hess->resize(ntets * 4);
//     }

//     tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
//         for (int i = range.begin(); i != range.end(); ++i) {
//             for (int j = 0; j < 4; j++) {
//                 tensors[i * 4 + j] = ComputeBasisProjectedTensorsPerTetFace(
//                     tet_facet_full_basis_.at(i).at(j), frames, is_kruskal, i, order,
//                     deriv ? &(*deriv)[i * 4 + j] : nullptr, hess ? &(*hess)[i * 4 + j] : nullptr);
//             }
//         }
//     });

//     // for (int i = 0; i != ntets; ++i) {
//     // 	for(int j = 0; j < 4; j++)
//     // 	{
//     // 		tensors[i*4 + j] = ComputeBasisProjectedTensorsPerTetFace(tet_facet_basis_.at(i).at(j), frames,
//     // is_kruskal, i, order, deriv ? &(*deriv)[i*4 + j] : nullptr, hess ? &(*hess)[i*4 + j] : nullptr);
//     // 	}
//     // }

//     return tensors;
// }

// Compute the frames projected to each tet face and then lifted as krushkal tensors
std::vector<Eigen::VectorXd> MiNTEnergy::ComputeDualEdgeProjectedTensors(
    const Eigen::MatrixXd& frames, bool is_kruskal, int order, std::vector<Eigen::MatrixXd>* deriv,
    std::vector<std::vector<Eigen::MatrixXd>>* hess) {
    int ntets = mesh_->nTets();

    std::vector<Eigen::VectorXd> tensors(ntets * 4);
    if (deriv) {
        deriv->resize(ntets * 4);
    }
    if (hess) {
        hess->resize(ntets * 4);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            for (int j = 0; j < 4; j++) {
                tensors[i * 4 + j] = ComputeBasisProjectedTensorsPerTetFace(
                    tet_facet_dual_basis_.at(i).at(j), frames, is_kruskal, i, order,
                    deriv ? &(*deriv)[i * 4 + j] : nullptr, hess ? &(*hess)[i * 4 + j] : nullptr);
            }
        }
    });

    return tensors;
}

// Compute the Kruskal tensor of a frame field projected to a tet face
Eigen::VectorXd MiNTEnergy::ComputeBasisProjectedTensorsPerTetFace(const Eigen::MatrixXd& facet_basis,
                                                                   const Eigen::MatrixXd& frames, bool is_kruskal,
                                                                   int tet_id, int order, Eigen::MatrixXd* deriv,
                                                                   std::vector<Eigen::MatrixXd>* hess) {
    int dim = 3;
    assert(frames.cols() % 3 == 0);
    int frame_per_tet = frames.cols() / dim;
    assert(frame_per_tet == 3);
    int proj_dim = facet_basis.rows();

    int coeff_size = SymmetricTensor::GetCoeffSize(proj_dim, order);

    Eigen::VectorXd tensor = Eigen::VectorXd::Zero(coeff_size);
    if (deriv) {
        deriv->resize(coeff_size, frames.cols());
        deriv->setZero();
    }

    if (hess) {
        hess->resize(coeff_size, Eigen::MatrixXd::Zero(frames.cols(), frames.cols()));
    }

    for (int i = 0; i < frame_per_tet; i++) {
        Eigen::VectorXd vec = frames.row(tet_id).segment(dim * i, dim).transpose();
        Eigen::VectorXd vec_proj = facet_basis * vec;

        Eigen::MatrixXd deriv_i;
        std::vector<Eigen::MatrixXd> hess_i;
        Eigen::VectorXd tensor_i;

        if (is_kruskal) {
            if (vec_proj.norm() < 1e-12) {
                vec_proj += 1e-12 * Eigen::VectorXd::Random(proj_dim);
            }
            tensor_i = SymmetricTensor::GetKruskalCoefficients(vec_proj, order, deriv ? &deriv_i : nullptr,
                                                               hess ? &hess_i : nullptr);
        } else {
            tensor_i = SymmetricTensor::GetMomentCoefficients(vec_proj, order, deriv ? &deriv_i : nullptr,
                                                              hess ? &hess_i : nullptr);
        }
        tensor += tensor_i;
        if (deriv) {
            // Eigen::MatrixXd blah = deriv_i * facet_basis;
            deriv->block(0, dim * i, coeff_size, dim) += deriv_i * facet_basis;
        }
        if (hess) {
            for (int j = 0; j < hess_i.size(); j++) {
                (*hess)[j].block(dim * i, dim * i, dim, dim) += facet_basis.transpose() * hess_i[j] * facet_basis;
            }
            // std::cout << "tensor: " << tensor.transpose() << std::endl;
        }
    }

    return tensor;
}

// Compute the fradeco tensors
std::vector<Eigen::VectorXd> MiNTEnergy::ComputeNormalizedMomentTensors(
    const Eigen::MatrixXd& frames, int order, std::vector<Eigen::MatrixXd>* deriv,
    std::vector<std::vector<Eigen::MatrixXd>>* hess) {
    int ntets = mesh_->nTets();

    std::vector<Eigen::VectorXd> tensors(ntets);
    if (deriv) {
        deriv->resize(ntets);
    }
    if (hess) {
        hess->resize(ntets);
    }

    // tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
    //     for (int i = range.begin(); i != range.end(); ++i) {
    //         tensors[i] = ComputeNormalizedMomentTensorsPerTet(frames, i, order, deriv ? &(*deriv)[i] : nullptr,
    //                                                           hess ? &(*hess)[i] : nullptr);
    //     }
    // });

    for (int i = 0; i < ntets; i++) {
        tensors[i] = ComputeNormalizedMomentTensorsPerTet(frames, i, order, deriv ? &(*deriv)[i] : nullptr,
                                                          hess ? &(*hess)[i] : nullptr);
    }

    return tensors;
}

// Compute the fradeco tensor of a frame field on a tetrahedron
Eigen::VectorXd MiNTEnergy::ComputeNormalizedMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                                 Eigen::MatrixXd* deriv,
                                                                 std::vector<Eigen::MatrixXd>* hess) {
    int dim = 3;
    assert(frames.cols() % 3 == 0);
    int frame_per_tet = frames.cols() / dim;
    assert(frame_per_tet == 3);

    int coeff_size = SymmetricTensor::GetCoeffSize(dim, order);

    Eigen::VectorXd kruskal_tensor = Eigen::VectorXd::Zero(coeff_size);
    if (deriv) {
        deriv->resize(coeff_size, frames.cols());
        deriv->setZero();
    }

    if (hess) {
        hess->resize(coeff_size, Eigen::MatrixXd::Zero(frames.cols(), frames.cols()));
    }

    for (int i = 0; i < frame_per_tet; i++) {
        Eigen::VectorXd vec = frames.row(tet_id).segment(dim * i, dim).transpose();
        Eigen::MatrixXd deriv_i;
        std::vector<Eigen::MatrixXd> hess_i;
        Eigen::VectorXd kruskal_tensor_i = SymmetricTensor::GetNormalizedMomentCoefficients(
            vec, order, deriv ? &deriv_i : nullptr, hess ? &hess_i : nullptr);

        kruskal_tensor += kruskal_tensor_i;
        if (deriv) {
            deriv->block(0, dim * i, coeff_size, dim) += deriv_i;
        }
        if (hess) {
            for (int j = 0; j < hess_i.size(); j++) {
                (*hess)[j].block(dim * i, dim * i, dim, dim) += hess_i[j];
            }
        }
    }

    return kruskal_tensor;
}

///////////// Assemble Hessian /////////////
// Assemble form per tet computation
double MiNTEnergy::AssembleEnergyFromTets(
    const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv, std::vector<Eigen::Triplet<double>>* hess_triplets,
    bool is_PSD_proj,
    std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    int ntets = mesh_->nTets();

    // this is used for parallel computation
    std::vector<double> energy_per_tet(ntets);
    std::vector<Eigen::VectorXd> deriv_per_tet(ntets);
    std::vector<Eigen::MatrixXd> hess_per_tet(ntets);

    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            energy_per_tet[i] += energy_func(frames, i, deriv ? &(deriv_per_tet[i]) : nullptr,
                                             hess_triplets ? &(hess_per_tet[i]) : nullptr, is_PSD_proj);
        }
    });

    // for (int i = 0; i != ntets; ++i) {
    //     energy_per_tet[i] += energy_func(frames, i, deriv ? &(deriv_per_tet[i]) : nullptr,
    //                                      hess_triplets ? &(hess_per_tet[i]) : nullptr, is_PSD_proj);
    // }

    // if (show_energy_) {
    //     per_tet_energy_density_[i] += energy_per_tet[i];
    // }

    // assembly the energy, derivative and hessian
    double energy = 0;
    int nvars_per_tet = frames.cols();
    assert(nvars_per_tet == 9);

    for (int i = 0; i < ntets; i++) {
        energy += energy_per_tet[i];
        if (show_energy_) {
            per_tet_energy_density_[i] += energy_per_tet[i];
        }

        if (deriv) {
            deriv->segment<9>(i * nvars_per_tet) += deriv_per_tet[i].segment<9>(0);
        }

        if (hess_triplets) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 9; k++) {
                    hess_triplets->push_back(
                        Eigen::Triplet<double>(i * nvars_per_tet + j, i * nvars_per_tet + k, hess_per_tet[i](j, k)));
                }
            }
        }
    }

    return energy;
}

// Assemble from per face computation
double MiNTEnergy::AssembleEnergyFromFaces(
    const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv, std::vector<Eigen::Triplet<double>>* hess_triplets,
    bool is_PSD_proj,
    std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func) {
    int nfaces = mesh_->nFaces();
    std::vector<double> energy_per_face(nfaces);
    std::vector<Eigen::VectorXd> deriv_per_face(nfaces);
    std::vector<Eigen::MatrixXd> hess_per_face(nfaces);

    tbb::parallel_for(tbb::blocked_range<int>(0u, nfaces), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            energy_per_face[i] += energy_func(frames, i, deriv ? &(deriv_per_face[i]) : nullptr,
                                              hess_triplets ? &(hess_per_face[i]) : nullptr, is_PSD_proj);
        }
    });

    // for (int i = 0; i != nfaces; ++i) {
    //     energy_per_face[i] += energy_func(frames, i, deriv ? &(deriv_per_face[i]) : nullptr,
    //                                       hess_triplets ? &(hess_per_face[i]) : nullptr, is_PSD_proj);
    // }

    // assembly the energy, derivative and hessian
    double energy = 0;
    int nvars_per_tet = frames.cols();
    assert(nvars_per_tet == 9);

    for (int i = 0; i < nfaces; i++) {
        int tet_id0 = mesh_->faceTet(i, 0);
        int tet_id1 = mesh_->faceTet(i, 1);

        if (tet_id0 == -1 || tet_id1 == -1) {
            continue;
        }
        energy += energy_per_face[i];

        if (show_energy_) {
            per_tet_energy_density_[tet_id0] += energy_per_face[i];
            per_tet_energy_density_[tet_id1] += energy_per_face[i];
        }

        if (std::isnan(energy_per_face[i])) {
            std::cout << "Energy is nan" << std::endl;
            std::cout << "Face id: " << i << std::endl;
            std::cout << "Tet id: 0" << tet_id0 << std::endl;
            std::cout << "tet vertex 0: " << mesh_->tetVertex(tet_id0, 0) << " "
                      << V_->row(mesh_->tetVertex(tet_id0, 0)) << std::endl;
            std::cout << "tet vertex 1: " << mesh_->tetVertex(tet_id0, 1) << " "
                      << V_->row(mesh_->tetVertex(tet_id0, 1)) << std::endl;
            std::cout << "tet vertex 2: " << mesh_->tetVertex(tet_id0, 2) << " "
                      << V_->row(mesh_->tetVertex(tet_id0, 2)) << std::endl;
            std::cout << "tet vertex 3: " << mesh_->tetVertex(tet_id0, 3) << " "
                      << V_->row(mesh_->tetVertex(tet_id0, 3)) << std::endl;

            std::cout << "Tet id: 1" << tet_id1 << std::endl;
            std::cout << "tet vertex 0: " << mesh_->tetVertex(tet_id1, 0) << " "
                      << V_->row(mesh_->tetVertex(tet_id1, 0)) << std::endl;
            std::cout << "tet vertex 1: " << mesh_->tetVertex(tet_id1, 1) << " "
                      << V_->row(mesh_->tetVertex(tet_id1, 1)) << std::endl;
            std::cout << "tet vertex 2: " << mesh_->tetVertex(tet_id1, 2) << " "
                      << V_->row(mesh_->tetVertex(tet_id1, 2)) << std::endl;
            std::cout << "tet vertex 3: " << mesh_->tetVertex(tet_id1, 3) << " "
                      << V_->row(mesh_->tetVertex(tet_id1, 3)) << std::endl;

            // << " " << tet_id1 << std::endl;
            std::cout << "Frame 1: " << frames.row(tet_id0) << std::endl;
            std::cout << "Frame 2: " << frames.row(tet_id1) << std::endl;
        }

        if (deriv) {
            deriv->segment<9>(tet_id0 * nvars_per_tet) += deriv_per_face[i].segment<9>(0);
            deriv->segment<9>(tet_id1 * nvars_per_tet) += deriv_per_face[i].segment<9>(9);
        }

        if (hess_triplets) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 9; k++) {
                    hess_triplets->push_back(Eigen::Triplet<double>(
                        tet_id0 * nvars_per_tet + j, tet_id0 * nvars_per_tet + k, hess_per_face[i](j, k)));

                    hess_triplets->push_back(Eigen::Triplet<double>(
                        tet_id0 * nvars_per_tet + j, tet_id1 * nvars_per_tet + k, hess_per_face[i](j, 9 + k)));

                    hess_triplets->push_back(Eigen::Triplet<double>(
                        tet_id1 * nvars_per_tet + j, tet_id0 * nvars_per_tet + k, hess_per_face[i](9 + j, k)));

                    hess_triplets->push_back(Eigen::Triplet<double>(
                        tet_id1 * nvars_per_tet + j, tet_id1 * nvars_per_tet + k, hess_per_face[i](9 + j, 9 + k)));
                }
            }
        }
    }

    if (std::isnan(energy)) {
        throw std::runtime_error("Energy is nan");
    }
    // assert(!std::isnan(energy));

    return energy;
}

///////////// Integrability Energy //////////////
// Compute the integrability energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols()
// == 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeIntegrabilityEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                          std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                          bool is_PSD_proj, bool is_kruskal, double weight) {
    double energy = 0.0;

    // sanity check for precomputation
    if (!mesh_ || !V_ || dual_edge_lengths_.size() != mesh_->nFaces() || dual_face_volumes_.size() != mesh_->nFaces() ||
        tensor2D_ord_weight_.size() != 3) {
        Precompute();
    }

    Eigen::MatrixXd frames_curr = frames;

    // numerical hack.  really should be addressed in the tensor derivative code
    if (is_kruskal) {
        frames_curr += 1e-10 * Eigen::MatrixXd::Random(frames.rows(), frames.cols());
    }

    int prev_start = 0;
    for (int ord = 2; ord <= 6; ord += 2) {
        // for (int ord = 2; ord <= 2; ord += 2) {
        Eigen::VectorXd deriv_ord;
        // std::cout << " trip size " << hess_triplets->size() << std::endl;
        // std::cout << "order: " << ord << std::endl;

        if (hess_triplets) {
            prev_start = hess_triplets->size();
            // std::cout << "prev_start " << prev_start << std::endl;
        }

        energy += ComputeIntegrabilityEnergyByOrderWithTriplets(frames_curr, ord, deriv ? &deriv_ord : nullptr,
                                                                hess_triplets, is_PSD_proj, is_kruskal, weight);

        if (hess_triplets && false) {
            std::vector<Eigen::Triplet<double>> newHess(hess_triplets->begin() + prev_start, hess_triplets->end());

            Eigen::SparseMatrix<double> hess_tmp;   // = new Eigen::SparseMatrix<double>;
            hess_tmp.resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
            hess_tmp.setFromTriplets(newHess.begin(), newHess.end());

            Eigen::MatrixXd hess_dense = Eigen::MatrixXd(hess_tmp);

            std::string filename = "tmp_" + std::to_string(ord) + ".csv";
            std::ofstream file(filename);
            if (file.is_open()) {
                file << hess_dense;
                file.close();
            }
        }

        if (deriv) {
            if (ord == 2) {
                *deriv = deriv_ord;
            } else {
                *deriv += deriv_ord;
            }
        }
    }

    return energy;
}

// Compute the integrability energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols()
// == 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeIntegrabilityEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                              Eigen::SparseMatrix<double>* hess, bool is_PSD_proj, bool is_kruskal) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeIntegrabilityEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj, is_kruskal);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());

        // if (false)
        // {
        //     Eigen::SparseMatrix<double> hess_tmp;
        //     hess_tmp.resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        //     hess_tmp.setFromTriplets(T.begin(), T.end());
        //     Eigen::MatrixXd hess_dense = Eigen::MatrixXd(hess_tmp);

        //     std::string filename = "tmp_full.csv";
        //     std::ofstream file(filename);
        //     if (file.is_open()) {
        //         file << hess_dense;
        //         file.close();
        //     }
        //     // throw std::runtime_error("testing");

        // }
    }

    return energy;
}

// This is the L1 integrability penalty.

double MiNTEnergy::ComputeIntegrabilityEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                                     std::vector<Eigen::MatrixXd>* tensor_deriv,
                                                     std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                                     const Eigen::VectorXd& coeff_weights, int face_id,
                                                     Eigen::VectorXd* deriv, Eigen::MatrixXd* hess, bool is_PSD_proj) {
    int nvars = 18;   // 3 frames per tet, 3 dimensions per frame, and we have 2 adjacent tets
    if (deriv) {
        deriv->setZero(nvars);
        if (!tensor_deriv) {
            std::cout << "Warning: The tensor derivative is not provided" << std::endl;
        }
    }
    if (hess) {
        if (!tensor_deriv || !tensor_hess) {
            std::cout << "Warning: The tensor derivative or hessian is not provided" << std::endl;
        }
    }

    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0.;   // bnd tet
    }

    int ntets_orig = mesh_not_extended_->nTets();

    // Enforce across boundary faces if below is comment out
    // if (tet_id0 >= ntets_orig || tet_id1 >= ntets_orig) {
    //     // std::cout << ntets_orig << " " << tet_id0 << " " << tet_id1 << std::endl;

    //     return 0.;   // bnd tet
    // }

    // TODO: update this to use dual face area / dual edge length instead.
    // This is almost the same, but that one is more principled.
    // also this should be moved to the precompute step.
    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    // weight *= this->weights_.w_mint;


    int tetfaceidx0 = tet_id0 * 4 + mesh_->faceTetIndex(face_id, 0);
    int tetfaceidx1 = tet_id1 * 4 + mesh_->faceTetIndex(face_id, 1);

    Eigen::VectorXd tensor_diff = (tensors[tetfaceidx0] - tensors[tetfaceidx1]);
    // std::cout << "tensor_diff: " << tensor_diff.transpose() << " | id: " << tetfaceidx0 << " | pt: " <<  " "<<
    // tensors[tetfaceidx0].transpose() << " | id: " << tetfaceidx1 << " | pt: " << tensors[tetfaceidx1].transpose() <<
    // std::endl;

    double energy = 0;
    Eigen::MatrixXd top_left_block = Eigen::MatrixXd::Zero(9, 9);
    Eigen::MatrixXd bottom_right_block = Eigen::MatrixXd::Zero(9, 9);

    // L1 term

    double l1 = .0001 * this->weights_.w_smooth * 0.  ;
    double l2 = 1.0 * this->weights_.w_mint;

    for (int i = 0; i < tensor_diff.size(); i++) {
        double sign = (tensor_diff[i] == 0) ? 0 : (tensor_diff[i] > 0 ? 1 : -1); // Subgradient

        energy += coeff_weights[i] * sign * tensor_diff[i] * l1;

        if ((deriv || hess) && tensor_deriv) {
            Eigen::VectorXd tensor_diff_i_deriv(18);
            tensor_diff_i_deriv.setZero();
            tensor_diff_i_deriv.segment<9>(0) = (*tensor_deriv)[tetfaceidx0].row(i).transpose() * sign;
            tensor_diff_i_deriv.segment<9>(9) = -(*tensor_deriv)[tetfaceidx1].row(i).transpose() * sign;

            if (deriv) {
                (*deriv) += coeff_weights[i] * tensor_diff_i_deriv * l1;
            }

            if (hess && tensor_hess) {
                // if (i == 0) {
                //     (*hess) = coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                // } else {
                //     (*hess) += coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                // }

                (*hess) = coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose() * 0.;


                top_left_block += sign * coeff_weights[i] * (*tensor_hess)[tetfaceidx0][i] * l1;
                bottom_right_block -= sign * coeff_weights[i] * (*tensor_hess)[tetfaceidx1][i] * l1;
            }
        }
    }

    // L2 term
    for (int i = 0; i < tensor_diff.size(); i++) {
        energy += coeff_weights[i] * tensor_diff[i] * tensor_diff[i] * l2;

        if ((deriv || hess) && tensor_deriv) {
            Eigen::VectorXd tensor_diff_i_deriv(18);
            tensor_diff_i_deriv.setZero();
            tensor_diff_i_deriv.segment<9>(0) = (*tensor_deriv)[tetfaceidx0].row(i).transpose();
            tensor_diff_i_deriv.segment<9>(9) = -(*tensor_deriv)[tetfaceidx1].row(i).transpose();

            if (deriv) {
                (*deriv) += 2 * tensor_diff[i] * coeff_weights[i] * tensor_diff_i_deriv * l2;
            }

            if (hess && tensor_hess) {
                if (i == 0) {
                    (*hess) += 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose() * l2;
                } else {
                    (*hess) += 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose()  * l2;
                }

                top_left_block += 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tetfaceidx0][i] * l2;
                bottom_right_block -= 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tetfaceidx1][i] * l2;
            }
        }
    }

    energy *= weight;

    if (deriv) {
        (*deriv) *= weight;
    }

    if (hess) {
        if (is_PSD_proj) {
            //(*hess) = SelfAdjSPDProjection(*hess);
            hess->block<9, 9>(0, 0) += SelfAdjSPDProjection(top_left_block);
            hess->block<9, 9>(9, 9) += SelfAdjSPDProjection(bottom_right_block);
        } else {
            hess->block<9, 9>(0, 0) += top_left_block;
            hess->block<9, 9>(9, 9) += bottom_right_block;
        }
        (*hess) *= weight;

        // std::cout << "tensor weights: " << coeff_weights.transpose() << " weight: " << weight << " fid: " << face_id
        // << " t0: " << tet_id0 << " t1: " << tet_id1 << " e: " << energy << std::endl;
    }

    return energy;
}




// This is the L2 integrability penalty.
// /*

// Compute the integrability energy of the frame field per face by order
// this is almost the same as the smoothness energy except should use different weights and per face indexing
double MiNTEnergy::ComputeIntegrabilityEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                                     std::vector<Eigen::MatrixXd>* tensor_deriv,
                                                     std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                                     const Eigen::VectorXd& coeff_weights, int face_id,
                                                     Eigen::VectorXd* deriv, Eigen::MatrixXd* hess, bool is_PSD_proj) {
    int nvars = 18;   // 3 frames per tet, 3 dimensions per frame, and we have 2 adjacent tets
    if (deriv) {
        deriv->setZero(nvars);
        if (!tensor_deriv) {
            std::cout << "Warning: The tensor derivative is not provided" << std::endl;
        }
    }
    if (hess) {
        if (!tensor_deriv || !tensor_hess) {
            std::cout << "Warning: The tensor derivative or hessian is not provided" << std::endl;
        }
    }

    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0.;   // bnd tet
    }

    int ntets_orig = mesh_not_extended_->nTets();

    // if (tet_id0 >= ntets_orig || tet_id1 >= ntets_orig) {
    //     // std::cout << ntets_orig << " " << tet_id0 << " " << tet_id1 << std::endl;

    //     return 0.;   // bnd tet
    // }

    // TODO: update this to use dual face area / dual edge length instead.
    // This is almost the same, but that one is more principled.
    // also this should be moved to the precompute step.
    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    // weight *= this->weights_.w_mint;

    // if exit early there are issues with hessian assembly.
    // if (tet_id0 >= ntets_orig || tet_id1 >= ntets_orig) {
    //     weight = 0;
    // }

    int tetfaceidx0 = tet_id0 * 4 + mesh_->faceTetIndex(face_id, 0);
    int tetfaceidx1 = tet_id1 * 4 + mesh_->faceTetIndex(face_id, 1);

    Eigen::VectorXd tensor_diff = (tensors[tetfaceidx0] - tensors[tetfaceidx1]);
    // std::cout << "tensor_diff: " << tensor_diff.transpose() << " | id: " << tetfaceidx0 << " | pt: " <<  " "<<
    // tensors[tetfaceidx0].transpose() << " | id: " << tetfaceidx1 << " | pt: " << tensors[tetfaceidx1].transpose() <<
    // std::endl;

    double energy = 0;
    Eigen::MatrixXd top_left_block = Eigen::MatrixXd::Zero(9, 9);
    Eigen::MatrixXd bottom_right_block = Eigen::MatrixXd::Zero(9, 9);

    for (int i = 0; i < tensor_diff.size(); i++) {
        energy += coeff_weights[i] * tensor_diff[i] * tensor_diff[i];

        if ((deriv || hess) && tensor_deriv) {
            Eigen::VectorXd tensor_diff_i_deriv(18);
            tensor_diff_i_deriv.setZero();
            tensor_diff_i_deriv.segment<9>(0) = (*tensor_deriv)[tetfaceidx0].row(i).transpose();
            tensor_diff_i_deriv.segment<9>(9) = -(*tensor_deriv)[tetfaceidx1].row(i).transpose();

            if (deriv) {
                (*deriv) += 2 * tensor_diff[i] * coeff_weights[i] * tensor_diff_i_deriv;
            }

            if (hess && tensor_hess) {
                if (i == 0) {
                    (*hess) = 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                } else {
                    (*hess) += 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                }

                top_left_block += 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tetfaceidx0][i];
                bottom_right_block -= 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tetfaceidx1][i];
            }
        }
    }

    energy *= weight;

    if (deriv) {
        (*deriv) *= weight;
    }

    if (hess) {
        if (is_PSD_proj) {
            //(*hess) = SelfAdjSPDProjection(*hess);
            hess->block<9, 9>(0, 0) += SelfAdjSPDProjection(top_left_block);
            hess->block<9, 9>(9, 9) += SelfAdjSPDProjection(bottom_right_block);
        } else {
            hess->block<9, 9>(0, 0) += top_left_block;
            hess->block<9, 9>(9, 9) += bottom_right_block;
        }
        (*hess) *= weight;

        // std::cout << "tensor weights: " << coeff_weights.transpose() << " weight: " << weight << " fid: " << face_id
        // << " t0: " << tet_id0 << " t1: " << tet_id1 << " e: " << energy << std::endl;
    }

    return energy;
}

//

// END L2 integrability penalty.
/////////////////////////////////////////////////////////////////////

// Compute the primal integrability energy of the frame field per face by order.
// This is the thing we actually want to be zero.
double MiNTEnergy::ComputePrimalIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id) {
    int nvars = 18;   // 3 frames per tet, 3 dimensions per frame, and we have 2 adjacent tets

    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0;   // bnd tet
    }

    // TODO: update this to use dual face area / dual edge length instead.
    // This is almost the same, but that one is more principled.
    // also this should be moved to the precompute step.
    // also this shouldn't really matter because it should uniformly go to numerical zero.
    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    // weight *= this->weights_.w_mint;

    int tetfaceidx0 = tet_id0 * 4 + mesh_->faceTetIndex(face_id, 0);
    int tetfaceidx1 = tet_id1 * 4 + mesh_->faceTetIndex(face_id, 1);

    Eigen::MatrixXd face_basis0 = tet_facet_basis_.at(tet_id0).at(mesh_->faceTetIndex(face_id, 0));
    // Eigen::MatrixXd face_basis1 = tet_facet_basis_.at(tet_id1).at(mesh_->faceTetIndex(face_id,1));

    Eigen::Vector3d int_error;

    for (int i = 0; i < 3; i++) {
        Eigen::VectorXd vec0 = frames.row(tet_id0).segment(3 * i, 3);
        Eigen::VectorXd proj_vec0 = face_basis0 * vec0;

        double smallest_diff = 1e10;

        for (int j = 0; j < 3; j++) {
            Eigen::VectorXd vec1 = frames.row(tet_id1).segment(3 * j, 3);
            Eigen::VectorXd proj_vec1 = face_basis0 * vec1;
            if (proj_vec0.dot(proj_vec1) < 0.) {
                proj_vec1 = -proj_vec1;
            }
            double cur_diff = (proj_vec0 - proj_vec1).squaredNorm();

            if (cur_diff < smallest_diff) {
                smallest_diff = cur_diff;
            }
        }
        int_error(i) = smallest_diff;
    }

    return int_error.squaredNorm() * weight;
}

// Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeIntegrabilityEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                                 Eigen::VectorXd* deriv,
                                                                 std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                                 bool is_PSD_proj, bool is_kruskal, double weight) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    if (weight == -1.) {
        // weight = this->weights_.w_mint;
        double scale = this->weights_.getGlobalScale();
        double lambda_penalty = this->weights_.lambda_penalty;

        double scaled_part = scale * lambda_penalty * this->weights_.w_int_sym;
        double fixed_part = scale * this->weights_.w_int_sym_fixed_weight;

        weight = scaled_part + fixed_part;
    }

    // std::cout << " integrability order " << ord << std::endl;

    std::vector<Eigen::MatrixXd> tensor_deriv;
    std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;
    // bool is_kruskal = false;   // add this to function api to do "augmented lagrangian" style opt
    // bool is_kruskal = true; // add this to function api to do "augmented lagrangian" style opt

    std::vector<Eigen::VectorXd> tensors =
        ComputeFaceProjectedTensors(frames, is_kruskal, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
                                    hess_triplets ? &tensor_hess : nullptr);

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        double energy = ComputeIntegrabilityEnergyPerFace(tensors, (face_deriv || face_hess) ? &tensor_deriv : nullptr,
                                                          face_hess ? &tensor_hess : nullptr, tensor2D_ord_weight_[ord],
                                                          face_id, face_deriv, face_hess, is_face_PSD_proj);
        if (face_deriv) *face_deriv *= weight;
        if (face_hess) *face_hess *= weight;
        energy *= weight;

        // if (face_deriv) *face_deriv *= this->weights_.w_smooth;
        // if (face_hess) *face_hess *= this->weights_.w_smooth;
        // energy *= this->weights_.w_smooth;
        return energy;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

//

// Compute the integrable energy of the frame field by order
double MiNTEnergy::ComputeIntegrabilityEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv,
                                                     Eigen::SparseMatrix<double>* hess, bool is_PSD_proj,
                                                     bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    std::vector<Eigen::Triplet<double>> T;

    double energy =
        ComputeIntegrabilityEnergyByOrderWithTriplets(frames, ord, deriv, hess ? &T : nullptr, is_PSD_proj, is_kruskal);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// this isn't quite right because it doesn't measure integrability error across boundary.
double MiNTEnergy::ComputePrimalIntegrabilityEnergy(const Eigen::MatrixXd& frames) {
    double energy = 0;

    int ntets = mesh_->nFaces();
    std::vector<double> energy_per_tet(ntets);
    tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            // if (mesh_->faceTet(i, 0) == -1 || mesh_->faceTet(i, 1) == -1) {
            // 	continue;
            // }
            energy_per_tet[i] = ComputePrimalIntegrabilityEnergyPerFace(frames, i);
        }
    });

    for (int i = 0; i < ntets; i++) {
        energy += energy_per_tet[i];
    }

    return energy;
}

///////////// Divergence Free Smoothness Energy //////////////
// Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron

// Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeDivergenceEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                              Eigen::VectorXd* deriv,
                                                              std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                              bool is_PSD_proj, bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    // std::cout << " integrability order " << ord << std::endl;

    std::vector<Eigen::MatrixXd> tensor_deriv;
    std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;
    // bool is_kruskal = false;   // add this to function api to do "augmented lagrangian" style opt
    // bool is_kruskal = true; // add this to function api to do "augmented lagrangian" style opt
    // is_kruskal = false;

    std::vector<Eigen::VectorXd> tensors =
        ComputeDualEdgeProjectedTensors(frames, is_kruskal, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
                                        hess_triplets ? &tensor_hess : nullptr);

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        Eigen::VectorXd weight = Eigen::VectorXd::Ones(1);
        double energy = ComputeIntegrabilityEnergyPerFace(tensors, (face_deriv || face_hess) ? &tensor_deriv : nullptr,
                                                          face_hess ? &tensor_hess : nullptr, weight, face_id,
                                                          face_deriv, face_hess, is_face_PSD_proj);
        // TODO make new variables
        if (face_deriv) *face_deriv *= this->weights_.w_smooth * 1e-1;
        if (face_hess) *face_hess *= this->weights_.w_smooth * 1e-1;
        return energy *= this->weights_.w_smooth * 1e-1;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

double MiNTEnergy::ComputeDivergenceEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                       std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                       bool is_PSD_proj, bool is_kruskal) {
    double energy = 0.0;

    // sanity check for precomputation
    if (!mesh_ || !V_ || dual_edge_lengths_.size() != mesh_->nFaces() || dual_face_volumes_.size() != mesh_->nFaces() ||
        tensor2D_ord_weight_.size() != 3) {
        Precompute();
    }

    Eigen::MatrixXd frames_curr = frames;

    // numerical hack.  really should be addressed in the tensor derivative code
    if (is_kruskal) {
        frames_curr += 1e-10 * Eigen::MatrixXd::Random(frames.rows(), frames.cols());
    }

    int prev_start = 0;
    for (int ord = 2; ord <= 6; ord += 2) {
        // for (int ord = 2; ord <= 2; ord += 2) {
        Eigen::VectorXd deriv_ord;
        // std::cout << " trip size " << hess_triplets->size() << std::endl;
        // std::cout << "order: " << ord << std::endl;

        if (hess_triplets) {
            prev_start = hess_triplets->size();
            // std::cout << "prev_start " << prev_start << std::endl;
        }

        energy += ComputeDivergenceEnergyByOrderWithTriplets(frames_curr, ord, deriv ? &deriv_ord : nullptr,
                                                             hess_triplets, is_PSD_proj, is_kruskal);

        if (hess_triplets && false) {
            std::vector<Eigen::Triplet<double>> newHess(hess_triplets->begin() + prev_start, hess_triplets->end());

            Eigen::SparseMatrix<double> hess_tmp;   // = new Eigen::SparseMatrix<double>;
            hess_tmp.resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
            hess_tmp.setFromTriplets(newHess.begin(), newHess.end());

            Eigen::MatrixXd hess_dense = Eigen::MatrixXd(hess_tmp);

            std::string filename = "tmp_" + std::to_string(ord) + ".csv";
            std::ofstream file(filename);
            if (file.is_open()) {
                file << hess_dense;
                file.close();
            }
        }

        if (deriv) {
            if (ord == 2) {
                *deriv = deriv_ord;
            } else {
                *deriv += deriv_ord;
            }
        }
    }

    return energy;
}

// Compute the integrability energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols()
// == 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeDivergenceEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj, bool is_kruskal) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeDivergenceEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj, is_kruskal);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());

        // if (false)
        // {
        //     Eigen::SparseMatrix<double> hess_tmp;
        //     hess_tmp.resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        //     hess_tmp.setFromTriplets(T.begin(), T.end());
        //     Eigen::MatrixXd hess_dense = Eigen::MatrixXd(hess_tmp);

        //     std::string filename = "tmp_full.csv";
        //     std::ofstream file(filename);
        //     if (file.is_open()) {
        //         file << hess_dense;
        //         file.close();
        //     }
        //     // throw std::runtime_error("testing");

        // }
    }

    return energy;
}

// // Compute the primal integrability energy of the frame field per face by order.
// // This is the thing we actually want to be zero.
// double MiNTEnergy::ComputePrimalIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id) {
//     int nvars = 18;   // 3 frames per tet, 3 dimensions per frame, and we have 2 adjacent tets

//     int tet_id0 = mesh_->faceTet(face_id, 0);
//     int tet_id1 = mesh_->faceTet(face_id, 1);
//     if (tet_id0 == -1 || tet_id1 == -1) {
//         return 0;   // bnd tet
//     }

//     // TODO: update this to use dual face area / dual edge length instead.
//     // This is almost the same, but that one is more principled.
//     // also this should be moved to the precompute step.
//     // also this shouldn't really matter because it should uniformly go to numerical zero.
//     double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
//     // weight *= this->weights_.w_mint;

//     int tetfaceidx0 = tet_id0 * 4 + mesh_->faceTetIndex(face_id, 0);
//     int tetfaceidx1 = tet_id1 * 4 + mesh_->faceTetIndex(face_id, 1);

//     Eigen::MatrixXd face_basis0 = tet_facet_basis_.at(tet_id0).at(mesh_->faceTetIndex(face_id, 0));
//     // Eigen::MatrixXd face_basis1 = tet_facet_basis_.at(tet_id1).at(mesh_->faceTetIndex(face_id,1));

//     Eigen::Vector3d int_error;

//     for (int i = 0; i < 3; i++) {
//         Eigen::VectorXd vec0 = frames.row(tet_id0).segment(3 * i, 3);
//         Eigen::VectorXd proj_vec0 = face_basis0 * vec0;

//         double smallest_diff = 1e10;

//         for (int j = 0; j < 3; j++) {
//             Eigen::VectorXd vec1 = frames.row(tet_id1).segment(3 * j, 3);
//             Eigen::VectorXd proj_vec1 = face_basis0 * vec1;
//             if (proj_vec0.dot(proj_vec1) < 0.) {
//                 proj_vec1 = -proj_vec1;
//             }
//             double cur_diff = (proj_vec0 - proj_vec1).squaredNorm();

//             if (cur_diff < smallest_diff) {
//                 smallest_diff = cur_diff;
//             }
//         }
//         int_error(i) = smallest_diff;
//     }

//     return int_error.squaredNorm() * weight;
// }

// Compute the integrable energy of the frame field by order
double MiNTEnergy::ComputeDivergenceEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv,
                                                  Eigen::SparseMatrix<double>* hess, bool is_PSD_proj,
                                                  bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    std::vector<Eigen::Triplet<double>> T;

    double energy =
        ComputeDivergenceEnergyByOrderWithTriplets(frames, ord, deriv, hess ? &T : nullptr, is_PSD_proj, is_kruskal);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// // this isn't quite right because it doesn't measure integrability error across boundary.
// double MiNTEnergy::ComputePrimalIntegrabilityEnergy(const Eigen::MatrixXd& frames) {
//     double energy = 0;

//     int ntets = mesh_->nFaces();
//     std::vector<double> energy_per_tet(ntets);
//     tbb::parallel_for(tbb::blocked_range<int>(0u, ntets), [&](const tbb::blocked_range<int>& range) {
//         for (int i = range.begin(); i != range.end(); ++i) {
//             // if (mesh_->faceTet(i, 0) == -1 || mesh_->faceTet(i, 1) == -1) {
//             // 	continue;
//             // }
//             energy_per_tet[i] = ComputePrimalIntegrabilityEnergyPerFace(frames, i);
//         }
//     });

//     for (int i = 0; i < ntets; i++) {
//         energy += energy_per_tet[i];
//     }

//     return energy;
// }

///////////// Hodge Laplacian Smoothness Energy //////////////
// Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeHodgeLaplacianSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames,
                                                                     Eigen::VectorXd* deriv,
                                                                     std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                                     bool is_PSD_proj) {
    double energy = 0.0;

    // sanity check for precomputation
    if (!mesh_ || !V_ || dual_edge_lengths_.size() != mesh_->nFaces() || dual_face_volumes_.size() != mesh_->nFaces() ||
        tensor_ord_weight_.size() != 3) {
        Precompute();
    }

    for (int ord = 2; ord <= 6; ord += 2) {
        Eigen::VectorXd deriv_ord;

        energy += ComputeHodgeLaplacianEnergyByOrderWithTriplets(frames, ord, deriv ? &deriv_ord : nullptr,
                                                                 hess_triplets, is_PSD_proj);

        if (deriv) {
            if (ord == 2) {
                *deriv = deriv_ord;
            } else {
                *deriv += deriv_ord;
            }
        }
    }

    return energy;
}

// Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeHodgeLaplacianSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                         Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeHodgeLaplacianSmoothnessEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// This term computes the energy || div ||^2  + || curl ||^2
// this can be shown to be equivilant to minimizing the hodge laplacian.
// In this function we implement the symmetric extension of the former formulation
// over "krushkal tensors".
// This is a very similar energy to the smoothness energy, except that one is effectively
// the smoothness term using the combinatorial laplacian
// whereas this implements smoothness relative to minimizing the hodge laplacian.

// The mathematical argument is as follows:
// The Hodge laplacian effectively measures smoothness in the sense that the only zeros
// of the vector version are harmonic functions.

// This property does not hold for the energy in "ComputeSmoothnessEnergy"

// TODO: write test cases that compare this with the smoothness energy.
double MiNTEnergy::ComputeHodgeLaplacianEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                                  Eigen::VectorXd* deriv,
                                                                  std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                                  bool is_PSD_proj) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    // std::cout << " integrability order " << ord << std::endl;

    // std::vector<Eigen::MatrixXd> tensor_deriv;
    // std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;

    std::vector<Eigen::MatrixXd> tensor_int_deriv;
    std::vector<Eigen::MatrixXd> tensor_div_deriv;
    std::vector<std::vector<Eigen::MatrixXd>> tensor_int_hess;
    std::vector<std::vector<Eigen::MatrixXd>> tensor_div_hess;

    // bool is_kruskal = false;   // add this to function api to do "augmented lagrangian" style opt
    bool is_kruskal = true;   // we formulate the hodge laplacian relative to the "kruskal" lifting map.

    std::vector<Eigen::VectorXd> tensors_int =
        ComputeFaceProjectedTensors(frames, is_kruskal, ord, (deriv || hess_triplets) ? &tensor_int_deriv : nullptr,
                                    hess_triplets ? &tensor_int_hess : nullptr);

    std::vector<Eigen::VectorXd> tensors_div =
        ComputeDualEdgeProjectedTensors(frames, is_kruskal, ord, (deriv || hess_triplets) ? &tensor_div_deriv : nullptr,
                                        hess_triplets ? &tensor_div_hess : nullptr);

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        Eigen::VectorXd face_int_deriv;
        Eigen::MatrixXd face_int_hess;

        Eigen::VectorXd face_div_deriv;
        Eigen::MatrixXd face_div_hess;

        double int_obj = ComputeIntegrabilityEnergyPerFace(
            tensors_int, (face_deriv || face_hess) ? &tensor_int_deriv : nullptr,
            face_hess ? &tensor_int_hess : nullptr, tensor2D_ord_weight_[ord], face_id,
            face_deriv ? &face_int_deriv : nullptr, face_hess ? &face_int_hess : nullptr, is_face_PSD_proj);

        Eigen::VectorXd weight = Eigen::VectorXd::Ones(1);
        double div_obj = ComputeIntegrabilityEnergyPerFace(
            tensors_div, (face_deriv || face_hess) ? &tensor_div_deriv : nullptr,
            face_hess ? &tensor_div_hess : nullptr, weight, face_id, face_deriv ? &face_div_deriv : nullptr,
            face_hess ? &face_div_hess : nullptr, is_face_PSD_proj);

        if (face_deriv) {
            *face_deriv = face_int_deriv + face_div_deriv;
        }

        if (face_hess) {
            *face_hess = face_int_hess + face_div_hess;
        }

        // double int_obj = ComputeIntegrabilityEnergyPerFace(tensors_int, (face_deriv || face_hess) ? &tensor_deriv :
        // nullptr,
        //                                          face_hess ? &tensor_hess : nullptr, tensor2D_ord_weight_[ord],
        //                                          face_id, face_deriv, face_hess, is_face_PSD_proj);

        return int_obj + div_obj;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

///////////// Combed Integrability Energy //////////////
double MiNTEnergy::ComputeCombedIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id,
                                                           Eigen::VectorXd* deriv, Eigen::MatrixXd* hess,
                                                           bool is_PSD_proj) {
    // Retrieve tetrahedron indices connected by the face
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0;   // Boundary tetrahedron, no energy contribution
    }

    int dim = 3;

    // Calculate weight based on mesh properties and a predefined weight parameter
    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    // weight *= this->weights_.w_mint;

    double scale = this->weights_.getGlobalScale();
    double lambda_penalty = this->weights_.lambda_penalty;

    double scaled_part = scale * lambda_penalty * this->weights_.w_int_combed;
    double fixed_part = scale * this->weights_.w_int_combed_fixed_weight;

    weight *= scaled_part + fixed_part;

    // Extract the frames corresponding to the two tetrahedra
    Eigen::VectorXd v = frames.row(tet_id0);
    Eigen::VectorXd u = frames.row(tet_id1);

    // Retrieve and expand the face basis
    Eigen::MatrixXd face_basis0 = tet_facet_basis_.at(tet_id0).at(mesh_->faceTetIndex(face_id, 0));
    // face_basis0 = Eigen::MatrixXd::Identity(3, 3);

    // Eigen::MatrixXd basis_expand = Eigen::kroneckerProduct(face_basis0, Eigen::MatrixXd::Identity(dim, dim));

    Eigen::MatrixXd basis_expand = Eigen::MatrixXd::Zero(6, 9);
    basis_expand.block<2, 3>(0, 0) = face_basis0;
    basis_expand.block<2, 3>(2, 3) = face_basis0;
    basis_expand.block<2, 3>(4, 6) = face_basis0;

    // Compute the optimal permutation matrix between the frames
    // Eigen::MatrixXd P = ComputeOptimalPermutationMatrix(v, u, face_basis0);
    Eigen::MatrixXd P = permutation_as_integrable_as_possible_[face_id];

    // Apply permutation to the second frame vector
    Eigen::VectorXd pu = P * u;

    // Compute the difference vector in the expanded basis
    Eigen::VectorXd diff = basis_expand * (v - pu);

    // Compute the energy as the squared norm of the difference vector
    double energy = diff.squaredNorm() * weight;

    // If derivative (gradient) is requested
    if (deriv) {
        deriv->setZero(18);

        // Compute the gradient with respect to the first frame (v)
        deriv->segment<9>(0) = 2.0 * weight * basis_expand.transpose() * diff;

        // Compute the gradient with respect to the second frame (u)
        deriv->segment<9>(9) = -2.0 * weight * P.transpose() * basis_expand.transpose() * diff;
    }

    // If Hessian is requested
    if (hess) {
        hess->setZero(18, 18);

        // Compute the Hessian with respect to the first frame (v)
        hess->block<9, 9>(0, 0) = 2.0 * weight * basis_expand.transpose() * basis_expand;

        // Compute the Hessian with respect to the second frame (u)
        hess->block<9, 9>(9, 9) = 2.0 * weight * P.transpose() * basis_expand.transpose() * basis_expand * P;

        // Compute the mixed Hessians (cross terms)
        hess->block<9, 9>(0, 9) = -2.0 * weight * basis_expand.transpose() * basis_expand * P;
        hess->block<9, 9>(9, 0) = hess->block<9, 9>(0, 9).transpose();   // Transpose for symmetry
    }

    // Optionally project the Hessian to be positive semi-definite (PSD)
    if (is_PSD_proj && hess) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*hess);
        Eigen::VectorXd eigvals = es.eigenvalues();
        Eigen::MatrixXd eigvecs = es.eigenvectors();

        // Project the eigenvalues to be non-negative
        Eigen::VectorXd eigvals_proj = eigvals.cwiseMax(0.0);

        // Reconstruct the Hessian from the projected eigenvalues
        *hess = eigvecs * eigvals_proj.asDiagonal() * eigvecs.transpose();
    }

    return energy;
}

// Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeCombedIntegrabilityEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                                std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                                bool is_PSD_proj, bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    // initialize permutation matricies across each edge.
    // SetPermutationToIdentity();

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        double energy =
            ComputeCombedIntegrabilityEnergyPerFace(frames, face_id, face_deriv, face_hess, is_face_PSD_proj);

        return energy;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

///////////// Combed Smoothness Energy //////////////

// ComputeCombedSmoothnessEnergyWithTriplets

double MiNTEnergy::ComputeCombedSmoothnessEnergyPerFace(const Eigen::MatrixXd& frames, int face_id,
                                                        Eigen::VectorXd* deriv, Eigen::MatrixXd* hess,
                                                        bool is_PSD_proj) {
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0;   // Boundary tetrahedron, no energy contribution
    }

    int dim = 3;

    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    double scale = this->weights_.getGlobalScale();
    weight *= (this->weights_.w_smooth_combed * scale + this->weights_.w_smooth_combed_precon * scale * scale);

    Eigen::VectorXd v = frames.row(tet_id0);
    Eigen::VectorXd u = frames.row(tet_id1);

    Eigen::MatrixXd face_basis0;

    Eigen::MatrixXd P;
    if (this->weights_.b_smooth_asap_combing) {
        P = permutation_as_smooth_as_possible_[face_id];
    } else {
        P = permutation_as_integrable_as_possible_[face_id];
    }

    // Build the permutation matrix P
    // Eigen::MatrixXd P = permutation_as_smooth_as_possible_[face_id];

    // Apply permutation to frame u
    Eigen::VectorXd pu = P * u;

    // Compute the differences
    Eigen::VectorXd diff = v - pu;

    // Compute energy
    double energy = diff.squaredNorm() * weight;

    if (deriv) {
        // Compute derivatives (gradient)
        deriv->setZero(18);

        // Derivative with respect to the first tetrahedron's frames (v)
        deriv->segment<9>(0) = 2.0 * weight * diff;

        // Derivative with respect to the second tetrahedron's frames (u)
        deriv->segment<9>(9) = -2.0 * weight * (P.transpose() * diff);
    }

    if (hess) {
        // Compute Hessian
        hess->setZero(18, 18);

        // Hessian with respect to the first tetrahedron's frames (v)
        hess->block<9, 9>(0, 0) = 2.0 * weight * Eigen::MatrixXd::Identity(9, 9);

        // Hessian with respect to the second tetrahedron's frames (u)
        hess->block<9, 9>(9, 9) = 2.0 * weight * P.transpose() * P;

        // Mixed derivatives (cross terms)
        hess->block<9, 9>(0, 9) = -2.0 * weight * P;
        hess->block<9, 9>(9, 0) = -2.0 * weight * P.transpose();
    }

    // Optionally project Hessian to be positive semi-definite (PSD)
    if (is_PSD_proj && hess) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*hess);
        Eigen::VectorXd eigvals = es.eigenvalues();
        Eigen::MatrixXd eigvecs = es.eigenvectors();

        // Project the eigenvalues to be non-negative
        Eigen::VectorXd eigvals_proj = eigvals.cwiseMax(0.0);

        // Reconstruct the Hessian
        *hess = eigvecs * eigvals_proj.asDiagonal() * eigvecs.transpose();
    }

    return energy;
}

// Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeCombedSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                             std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                             bool is_PSD_proj, bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        double energy = ComputeCombedSmoothnessEnergyPerFace(frames, face_id, face_deriv, face_hess, is_face_PSD_proj);

        return energy;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

///////////// Smoothness Energy //////////////
// Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                       std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                       bool is_PSD_proj, bool is_kruskal) {
    double energy = 0.0;

    // sanity check for precomputation
    if (!mesh_ || !V_ || dual_edge_lengths_.size() != mesh_->nFaces() || dual_face_volumes_.size() != mesh_->nFaces() ||
        tensor_ord_weight_.size() != 3) {
        Precompute();
    }

    for (int ord = 2; ord <= 6; ord += 2) {
        Eigen::VectorXd deriv_ord;
        energy += ComputeSmoothnessEnergyByOrderWithTriplets(frames, ord, deriv ? &deriv_ord : nullptr, hess_triplets,
                                                             is_PSD_proj, is_kruskal);

        if (deriv) {
            if (ord == 2) {
                *deriv = deriv_ord;
            } else {
                *deriv += deriv_ord;
            }
        }
    }

    return energy;
}

// Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeSmoothnessEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

double MiNTEnergy::ComputeCombedSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                 Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeCombedSmoothnessEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

double MiNTEnergy::ComputeCombedIntegrabilityEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                    Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeCombedIntegrabilityEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the smoothness energy of the frame field per face by order
double MiNTEnergy::ComputeSmoothnessEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                                  std::vector<Eigen::MatrixXd>* tensor_deriv,
                                                  std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                                  const Eigen::VectorXd& coeff_weights, int face_id,
                                                  Eigen::VectorXd* deriv, Eigen::MatrixXd* hess, bool is_PSD_proj) {
    int nvars = 18;   // 3 frames per tet, 3 dimensions per frame, and we have 2 adjacent tets
    if (deriv) {
        deriv->setZero(nvars);
        if (!tensor_deriv) {
            std::cout << "Warning: The tensor derivative is not provided" << std::endl;
        }
    }
    if (hess) {
        hess->setZero(18, 18);

        if (!tensor_deriv || !tensor_hess) {
            std::cout << "Warning: The tensor derivative or hessian is not provided" << std::endl;
        }
    }

    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);
    if (tet_id0 == -1 || tet_id1 == -1) {
        return 0;   // bnd tet
    }

    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);

    double scale = this->weights_.getGlobalScale();
    weight *= (this->weights_.w_smooth_sym * scale + this->weights_.w_smooth_sym_precon * scale * scale);

    // if (this->weights_.w_mint > .5) {
    //     weight = 1e-8;
    // }
    // weight = 1;

    Eigen::VectorXd tensor_diff = (tensors[tet_id0] - tensors[tet_id1]);

    double energy = 0;
    Eigen::MatrixXd top_left_block = Eigen::MatrixXd::Zero(9, 9);
    Eigen::MatrixXd bottom_right_block = Eigen::MatrixXd::Zero(9, 9);

    for (int i = 0; i < tensor_diff.size(); i++) {
        energy += coeff_weights[i] * tensor_diff[i] * tensor_diff[i];

        if ((deriv || hess) && tensor_deriv) {
            Eigen::VectorXd tensor_diff_i_deriv(18);
            tensor_diff_i_deriv.setZero();
            tensor_diff_i_deriv.segment<9>(0) = (*tensor_deriv)[tet_id0].row(i).transpose();
            tensor_diff_i_deriv.segment<9>(9) = -(*tensor_deriv)[tet_id1].row(i).transpose();

            if (deriv) {
                (*deriv) += 2 * tensor_diff[i] * coeff_weights[i] * tensor_diff_i_deriv;
            }

            if (hess && tensor_hess) {
                if (i == 0) {
                    (*hess) = 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                } else {
                    (*hess) += 2 * coeff_weights[i] * tensor_diff_i_deriv * tensor_diff_i_deriv.transpose();
                }

                top_left_block += 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tet_id0][i];
                bottom_right_block -= 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tet_id1][i];

                // hess->block<9, 9>(0, 0) += 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tet_id0][i];
                // hess->block<9, 9>(9, 9) -= 2 * tensor_diff[i] * coeff_weights[i] * (*tensor_hess)[tet_id1][i];
            }
        }
    }

    energy *= weight;

    if (deriv) {
        (*deriv) *= weight;
    }

    if (hess) {
        if (is_PSD_proj) {
            //(*hess) = SelfAdjSPDProjection(*hess);
            hess->block<9, 9>(0, 0) += SelfAdjSPDProjection(top_left_block);
            hess->block<9, 9>(9, 9) += SelfAdjSPDProjection(bottom_right_block);
        } else {
            hess->block<9, 9>(0, 0) += top_left_block;
            hess->block<9, 9>(9, 9) += bottom_right_block;
        }
        (*hess) *= weight;
    }

    ////
    // Little test of this function.
    ////

    // if ((tet_id0 == 321 || tet_id1 == 321) && (tet_id0 == 1844 || tet_id1 == 1844)) {
    //     std::cout << "tet_id0: " << tet_id0 << " tet_id1: " << tet_id1 << " face_id: " << face_id << std::endl;
    //     std::cout << coeff_weights.transpose() << std::endl;
    //     std::cout << tensors[tet_id0].transpose() << std::endl;
    //     std::cout << tensors[tet_id1].transpose() << std::endl;
    //     std::cout << "difference: " << tensor_diff.transpose() << std::endl;
    //     std::cout << "energy: " << energy << std::endl;
    //     // std::cout << "dual_face_volumes_: " << dual_face_volumes_[face_id] << " dual_edge_lengths_: " <<
    //     // dual_edge_lengths_[face_id] << std::endl;
    // }

    return energy;
}

// Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeSmoothnessEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                              Eigen::VectorXd* deriv,
                                                              std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                              bool is_PSD_proj, bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    std::vector<Eigen::MatrixXd> tensor_deriv;
    std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;
    std::vector<Eigen::VectorXd> tensors;

    if (is_kruskal) {
        tensors = ComputeKruskalTensors(frames, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
                                        hess_triplets ? &tensor_hess : nullptr);
    } else {
        tensors = ComputeMomentTensors(frames, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
                                       hess_triplets ? &tensor_hess : nullptr);
    }

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
                           Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
        double energy = ComputeSmoothnessEnergyPerFace(tensors, (face_deriv || face_hess) ? &tensor_deriv : nullptr,
                                                       face_hess ? &tensor_hess : nullptr, tensor_ord_weight_[ord],
                                                       face_id, face_deriv, face_hess, is_face_PSD_proj);

        return energy;
    };

    return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// // Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that
// we
// // always append to the triplets instead of clearing it.
// double MiNTEnergy::ComputeSmoothnessEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
//                                                               Eigen::VectorXd* deriv,
//                                                               std::vector<Eigen::Triplet<double>>* hess_triplets,
//                                                               bool is_PSD_proj, bool is_kruskal) {
//     if (deriv) {
//         deriv->setZero(frames.rows() * frames.cols());
//     }

//     std::vector<Eigen::MatrixXd> tensor_deriv;
//     std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;
//     std::vector<Eigen::VectorXd> tensors;

//     if (is_kruskal) {
//         tensors = ComputeKruskalTensors(frames, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
//                                         hess_triplets ? &tensor_hess : nullptr);
//     } else {
//         tensors = ComputeMomentTensors(frames, ord, (deriv || hess_triplets) ? &tensor_deriv : nullptr,
//                                        hess_triplets ? &tensor_hess : nullptr);
//     }

//     auto energy_func = [&](const Eigen::MatrixXd& input_frames, int face_id, Eigen::VectorXd* face_deriv,
//                            Eigen::MatrixXd* face_hess, bool is_face_PSD_proj) {
//         double energy = ComputeSmoothnessEnergyPerFace(tensors, (face_deriv || face_hess) ? &tensor_deriv : nullptr,
//                                                        face_hess ? &tensor_hess : nullptr, tensor_ord_weight_[ord],
//                                                        face_id, face_deriv, face_hess, is_face_PSD_proj);

//         return energy;
//     };

//     return AssembleEnergyFromFaces(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
// }

// Compute the smoothness energy of the frame field by order
double MiNTEnergy::ComputeSmoothnessEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv,
                                                  Eigen::SparseMatrix<double>* hess, bool is_PSD_proj,
                                                  bool is_kruskal) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    std::vector<Eigen::Triplet<double>> T;

    double energy =
        ComputeSmoothnessEnergyByOrderWithTriplets(frames, ord, deriv, hess ? &T : nullptr, is_PSD_proj, is_kruskal);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

///////////// Compute Deviation Energy //////////////
// Compute the deviation energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeDeviationEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                      std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                      bool is_PSD_proj) {
    double energy = 0.0;

    if (guided_tensor_coeffs_.empty()) {
        std::cout << "guided frames are not provided" << std::endl;

        if (deriv) {
            deriv->setZero(frames.rows() * frames.cols());
        }

        return energy;
    }

    // sanity check for precomputation
    if (!mesh_ || !V_ || dual_edge_lengths_.size() != mesh_->nFaces() || dual_face_volumes_.size() != mesh_->nFaces() ||
        tensor_ord_weight_.size() != 3) {
        Precompute();
    }

    for (int ord = 2; ord <= 6; ord += 2) {
        Eigen::VectorXd deriv_ord;
        energy += ComputeDeviationEnergyByOrderWithTriplets(frames, ord, deriv ? &deriv_ord : nullptr, hess_triplets,
                                                            is_PSD_proj);

        if (deriv) {
            if (ord == 2) {
                *deriv = deriv_ord;
            } else {
                *deriv += deriv_ord;
            }
        }
    }

    return energy;
}

// Compute the deviation energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols() ==
// 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeDeviationEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                          Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeDeviationEnergyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the deviation energy per tet by order
double MiNTEnergy::ComputeDeviationEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                Eigen::VectorXd* deriv, Eigen::MatrixXd* hess, bool is_PSD_proj) {
    Eigen::VectorXd tensor;
    Eigen::MatrixXd tensor_deriv;
    std::vector<Eigen::MatrixXd> tensor_hess;

    switch (this->weights_.fit_type) {
        case FitType::Moments:
            // tensor = SymmetricTensor::GetMomentCoefficients(frames.row(tet_id), order, deriv ? &tensor_deriv :
            // nullptr,
            //                                                 hess ? &tensor_hess : nullptr);

            tensor = ComputeMomentTensorsPerTet(frames, tet_id, order, (deriv || hess) ? &tensor_deriv : nullptr,
                                                hess ? &tensor_hess : nullptr);

            // std::cout << tensor << std::endl;

            break;
        case FitType::Direction:
            tensor = ComputeNormalizedMomentTensorsPerTet(
                frames, tet_id, order, (deriv || hess) ? &tensor_deriv : nullptr, hess ? &tensor_hess : nullptr);
            break;
        case FitType::Kruskal:

            tensor = ComputeKruskalTensorsPerTet(frames, tet_id, order, (deriv || hess) ? &tensor_deriv : nullptr,
                                                 hess ? &tensor_hess : nullptr);
            // SymmetricTensor::GetMomentCoefficients(frames.row(tet_id), order, deriv ? &tensor_deriv : nullptr,
            //                                                           hess ? &tensor_hess : nullptr);
            // throw std::runtime_error("Kruskal tensor is not implemented here yet");
            break;
        default:
            break;
    }

    Eigen::VectorXd& guided_tensor = guided_tensor_coeffs_[order][tet_id];
    Eigen::VectorXd& tensor_weight = tensor_ord_weight_[order];
    double scale = this->weights_.getGlobalScale();

    double base_weight = tet_volumes_[tet_id] / ave_edge_len_ / ave_edge_len_;   // make this unit of length
                                                                                 // weight *= this->weights_.w_fit;

    base_weight = base_weight * scale * (this->weights_.w_fit + scale * this->weights_.w_fit_precon);

    // if (tet_id > mesh_not_extended_->nTets()) {
    //     // base_weight *= 1e-1 * this->weights_.w_bound_rescale;
    //     // tensor = ComputeMomentTensorsPerTet(frames, tet_id, order, (deriv || hess) ? &tensor_deriv : nullptr,
    //     //                                     hess ? &tensor_hess : nullptr);
    //     // base_weight *= 1.;

    //     tensor = ComputeNormalizedMomentTensorsPerTet(frames, tet_id, order, (deriv || hess) ? &tensor_deriv :
    //     nullptr,
    //                                                   hess ? &tensor_hess : nullptr);
    //     base_weight *= .001;

    // } else {
    //     // base_weight *= 1e1 * 0;
    //     tensor = ComputeNormalizedMomentTensorsPerTet(frames, tet_id, order, (deriv || hess) ? &tensor_deriv :
    //     nullptr,
    //                                                   hess ? &tensor_hess : nullptr);
    //     base_weight *= .001 * 0;
    // }
    // base_weight *= .01;
    // base_weight = 0;

    double energy = 0;
    for (int i = 0; i < tensor.size(); i++) {
        double weight = base_weight * tensor_weight[i];

        energy += weight * (tensor[i] - guided_tensor[i]) * (tensor[i] - guided_tensor[i]);

        if (deriv) {
            if (i == 0) {
                (*deriv) = 2 * weight * (tensor[i] - guided_tensor[i]) * tensor_deriv.row(i).transpose();
            } else {
                (*deriv) += 2 * weight * (tensor[i] - guided_tensor[i]) * tensor_deriv.row(i).transpose();
            }
        }

        if (hess) {
            if (i == 0) {
                (*hess) = 2 * weight * (tensor[i] - guided_tensor[i]) * tensor_hess[i];
            } else {
                (*hess) += 2 * weight * (tensor[i] - guided_tensor[i]) * tensor_hess[i];
            }
            (*hess) += 2 * weight * tensor_deriv.row(i).transpose() * tensor_deriv.row(i);
        }
    }

    if (hess && is_PSD_proj) {
        (*hess) = SelfAdjSPDProjection(*hess);
    }

    return energy;
}

// Compute the deviation energy per tet by order
double MiNTEnergy::ComputeOrthogonalityEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                    Eigen::MatrixXd* hess, bool is_PSD_proj) {
    int dim = 3;
    Eigen::Vector3d v1 = frames.row(tet_id).segment(0, dim);
    Eigen::Vector3d v2 = frames.row(tet_id).segment(3, dim);
    Eigen::Vector3d v3 = frames.row(tet_id).segment(6, dim);

    double e12 = v1.dot(v2);
    double e13 = v1.dot(v3);
    double e23 = v2.dot(v3);

    if (e12 / (v1.norm() * v2.norm()) > .99) {
        v1 = v1 + Eigen::Vector3d::Random() * 1e-6;
    }

    if (e13 / (v1.norm() * v3.norm()) > .99) {
        v3 = v3 + Eigen::Vector3d::Random() * 1e-6;
    }

    if (e23 / (v2.norm() * v3.norm()) > .99) {
        v2 = v2 + Eigen::Vector3d::Random() * 1e-6;
    }

    double weight = tet_volumes_[tet_id] / ave_edge_len_ / ave_edge_len_;   // make this unit of length
                                                                            // weight *= this->weights_.w_fit;

    double scale = this->weights_.getGlobalScale();
    weight *= (this->weights_.w_orthog_constraint_weight + this->weights_.w_orthog * scale +
               this->weights_.w_orthog_precon * scale * scale);
    // if (tet_id > mesh_not_extended_->nTets()) {
    //     // weight *= 1e2 * this->weights_.w_bound_rescale * 0.;
    //     weight *= 1e-3 * this->weights_.w_bound_rescale;

    // } else {
    //     weight *= .1;
    // }
    // weight *= 10000.;

    double energy = 0.5 * (e12 * e12 + e13 * e13 + e23 * e23);
    energy *= weight;

    if (deriv) {
        deriv->setZero(frames.cols(), 1);

        deriv->block(0, 0, dim, 1) = e12 * v2 + e13 * v3;
        deriv->block(3, 0, dim, 1) = e12 * v1 + e23 * v3;
        deriv->block(6, 0, dim, 1) = e13 * v1 + e23 * v2;

        (*deriv) *= weight;
    }

    if (hess) {
        hess->setZero(frames.cols(), frames.cols());

        Eigen::Matrix3d H11 = v2 * v2.transpose() + v3 * v3.transpose();
        Eigen::Matrix3d H22 = v1 * v1.transpose() + v3 * v3.transpose();
        Eigen::Matrix3d H33 = v1 * v1.transpose() + v2 * v2.transpose();

        Eigen::Matrix3d H12 = e12 * Eigen::Matrix3d::Identity() + v2 * v1.transpose();
        Eigen::Matrix3d H13 = e13 * Eigen::Matrix3d::Identity() + v3 * v1.transpose();
        Eigen::Matrix3d H23 = e23 * Eigen::Matrix3d::Identity() + v3 * v2.transpose();

        hess->block(0, 0, dim, dim) = H11;
        hess->block(0, 3, dim, dim) = H12;
        hess->block(0, 6, dim, dim) = H13;
        hess->block(3, 0, dim, dim) = H12.transpose();
        hess->block(3, 3, dim, dim) = H22;
        hess->block(3, 6, dim, dim) = H23;
        hess->block(6, 0, dim, dim) = H13.transpose();
        hess->block(6, 3, dim, dim) = H23.transpose();
        hess->block(6, 6, dim, dim) = H33;

        (*hess) *= weight;
    }

    if (hess && is_PSD_proj) {
        (*hess) = SelfAdjSPDProjection(*hess);
    }

    return energy;
}

// Compute the deviation energy of the frame field by order.
double MiNTEnergy::ComputeDeviationEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv,
                                                 Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    std::vector<Eigen::Triplet<double>> T;

    double energy = ComputeDeviationEnergyByOrderWithTriplets(frames, ord, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the deviation energy of the frame field by order, but we return the hessian in triplet form. Notice that we
// always append to the triplets instead of clearing it.
double MiNTEnergy::ComputeDeviationEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                             Eigen::VectorXd* deriv,
                                                             std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                             bool is_PSD_proj) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        return ComputeDeviationEnergyPerTet(input_frames, tet_id, ord, tet_deriv, tet_hess, is_tet_PSD_proj);
    };
    return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

///////////// Unit Norm Energy //////////////
// Compute the unit norm penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
// frames per tetrahedron
double MiNTEnergy::ComputeUnitNormPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                      std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                      bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        double energy = ComputeUnitNormPenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
        return energy;
    };
    return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// Compute the unit norm penalty
double MiNTEnergy::ComputeUnitNormPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                          Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeUnitNormPenaltyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the unit norm penalty per tet
double MiNTEnergy::ComputeUnitNormPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id] / 3.0;     // divided by 3.0 for the averaging through there frames
    weight *= 1. / ave_edge_len_ / ave_edge_len_;   // make it unit of length, scaled by weight coeff

    double base_weight = weight;

    if (deriv) {
        deriv->resize(9);
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    double scale = this->weights_.getGlobalScale();

    for (int i = 0; i < 3; i++) {
        Eigen::Vector3d v = frames.row(tet_id).segment<3>(3 * i);
        double diff = v.squaredNorm() - 1;   // v.squaredNorm() computes v[0]*v[0] + v[1]*v[1] + v[2]*v[2]

        // This is to set different weights

        if (i == 0 && tet_id > mesh_not_extended_->nTets()) {
            weight = base_weight * (this->weights_.w_smooth * 10. * this->weights_.w_bound_rescale * 0.);   // + 1e-8
        } else if (i > 0 && tet_id > mesh_not_extended_->nTets()) {
            // weight = base_weight * this->weights_.w_smooth * 1e4 * this->weights_.w_bound_rescale;
            weight = base_weight * scale * (this->weights_.w_unit + scale * this->weights_.w_unit_precon);

        } else {
            // weight = base_weight * this->weights_.w_smooth * 1e4;
            weight = base_weight * scale * (this->weights_.w_unit + scale * this->weights_.w_unit_precon);
            // double mu = 1e-1;
            // weight = base_weight * this->weights_.w_smooth * mu / 2 * .5 * 0.;
        }

        energy += weight * diff * diff;

        if (deriv || hess) {
            Eigen::Vector3d grad_diff = 2 * v;   // Gradient of diff with respect to v

            if (deriv) {
                (*deriv).segment<3>(3 * i) =
                    2 * weight * diff * grad_diff;   // Gradient of the energy with respect to v
            }

            if (hess) {
                Eigen::Matrix3d hess_diff = 4 * Eigen::Matrix3d::Identity();   // Hessian of diff with respect to v

                // Hessian of the energy with respect to v
                Eigen::Matrix3d hess_energy = weight * (diff * hess_diff + 2 * grad_diff * grad_diff.transpose());

                hess->block<3, 3>(3 * i, 3 * i) = hess_energy;
            }
        }
    }

    if (hess && is_PSD_proj) {
        for (int i = 0; i < 3; i++) {
            hess->block<3, 3>(3 * i, 3 * i) = SelfAdjSPDProjection(hess->block<3, 3>(3 * i, 3 * i));
        }
    }
    return energy;
}

///////////// Unit Norm Boundary Energy //////////////
// Make the first vector of each boundary element align with a unit normal. We assume that frames.rows() ==
// mesh_->nTets() and frames.cols() == 9, that is 3 frames per tetrahedron
double MiNTEnergy::ComputeUnitNormalBoundaryPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                                std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                                bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        double energy =
            ComputeUnitNormalBoundaryPenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
        return energy;
    };

    if (deriv) {
        deriv->setZero(frames.rows() * frames.cols());
    }

    int ntets = mesh_not_extended_->nTets();
    int nbound = mesh_not_extended_->nBoundaryElements();

    if (frames.rows() < ntets + nbound) {
        std::cout << "Not enough boundary DOFS not enforcing anything" << frames.rows() << " vs " << ntets + nbound
                  << std::endl;
        return 0;
    }

    // this is used for parallel computation
    std::vector<double> energy_per_tet(nbound);
    std::vector<Eigen::VectorXd> deriv_per_tet(nbound);
    std::vector<Eigen::MatrixXd> hess_per_tet(nbound);

    tbb::parallel_for(tbb::blocked_range<int>(0u, nbound), [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
            energy_per_tet[i] += energy_func(frames, i + ntets, deriv ? &(deriv_per_tet[i]) : nullptr,
                                             hess_triplets ? &(hess_per_tet[i]) : nullptr, is_PSD_proj);
        }
    });

    // assembly the energy, derivative and hessian
    double energy = 0;
    int nvars_per_tet = frames.cols();
    assert(nvars_per_tet == 9);

    for (int i = 0; i < nbound; i++) {
        energy += energy_per_tet[i];

        if (deriv) {
            deriv->segment<9>((i + ntets) * nvars_per_tet) += deriv_per_tet[i].segment<9>(0);
        }

        if (hess_triplets) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 9; k++) {
                    hess_triplets->push_back(Eigen::Triplet<double>(
                        (i + ntets) * nvars_per_tet + j, (i + ntets) * nvars_per_tet + k, hess_per_tet[i](j, k)));
                }
            }
        }
    }

    return energy;
    // return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// Compute the unit norm penalty
double MiNTEnergy::ComputeUnitNormalBoundaryPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                    Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeUnitNormalBoundaryPenaltyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the unit norm penalty per tet
double MiNTEnergy::ComputeUnitNormalBoundaryPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                                          Eigen::VectorXd* deriv, Eigen::MatrixXd* hess,
                                                          bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;

    double weight = tet_volumes_[tet_id];
    weight *= ((this->weights_.w_bound) / ave_edge_len_ / ave_edge_len_);   // make this unit of length, scaled by coeff

    // weight *= (this->weights_.w_bound / ave_edge_len_ / ave_edge_len_);   // make this unit of length, scaled by
    // coeff

    // This rescales the boundary weight so that when fields are unit, the area and volume terms in the objective scale
    // uniformly under refinement.
    weight *= this->weights_.w_bound_rescale * 0.;
    // weight *= (10. / ave_edge_len_ / ave_edge_len_);

    double orthog_weight = tet_volumes_[tet_id] * (10. / ave_edge_len_ / ave_edge_len_);
    // orthog_weight *= this->weights_.w_bound_rescale;
    // orthog_weight *= this->weights_.w_smooth;
    //  double orthog_weight = weight * 1000;

    if (deriv) {
        deriv->resize(9);
        deriv->setZero();
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    if (tet_id < mesh_not_extended_->nTets()) {
        return 0;
    }

    for (int i = 0; i < 1; i++) {
        Eigen::Vector3d v = frames.row(tet_id).segment<3>(3 * i);
        Eigen::Vector3d targ = boundary_normals_.row(tet_id - mesh_not_extended_->nTets());
        double diff = (v - targ).dot(v - targ);
        energy += weight * diff;

        Eigen::Vector3d b1 = boundary_b1_.row(tet_id - mesh_not_extended_->nTets());
        Eigen::Vector3d b2 = boundary_b2_.row(tet_id - mesh_not_extended_->nTets());

        energy += orthog_weight * (v.dot(b1) * v.dot(b1) + v.dot(b2) * v.dot(b2));

        if (deriv || hess) {
            Eigen::Vector3d grad_diff = 2 * (v - targ);
            Eigen::Vector3d grad_orthog = 2 * orthog_weight * (v.dot(b1) * b1 + v.dot(b2) * b2);

            if (deriv) {
                (*deriv).segment<3>(3 * i) += grad_diff * weight;
                (*deriv).segment<3>(3 * i) += grad_orthog;
            }

            if (hess) {
                Eigen::Matrix3d hess_diff = 2 * weight * Eigen::Matrix3d::Identity();
                Eigen::Matrix3d hess_orthog = 2 * orthog_weight * (b1 * b1.transpose() + b2 * b2.transpose());

                // Add diagonal terms for diff part
                hess->block<3, 3>(3 * i, 3 * i) += hess_diff;

                // Add full matrix for orthog part
                hess->block<3, 3>(3 * i, 3 * i) += hess_orthog;
            }
        }
    }

    if (hess && is_PSD_proj) {
        for (int i = 0; i < 1; i++) {
            hess->block<3, 3>(3 * i, 3 * i) = SelfAdjSPDProjection(hess->block<3, 3>(3 * i, 3 * i));
        }
    }
    return energy;
}

///////////// Unit Determinant Energy //////////////
// Compute the unit determinant penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is
// 3 frames per tetrahedron
double MiNTEnergy::ComputeDeterminantPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                         std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                         bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        return ComputeDeterminantPenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
    };
    return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// Compute the scaled jacobian penalty
double MiNTEnergy::ComputeDeterminantPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                             Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeDeterminantPenaltyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

///////////// Scaled Jacobian Energy //////////////
// Compute the scaled jacobian penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
// frames per tetrahedron
double MiNTEnergy::ComputeScaledJacobianPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                            std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                            bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        // return ComputeDeviationEnergyPerTet(input_frames, tet_id, 2, tet_deriv, tet_hess, is_tet_PSD_proj);
        return ComputeOrthogonalityEnergyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);

        // return ComputeScaledJacobianPenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
        // return ComputeScaledJacobianInversePenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
    };
    return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// Compute the scaled jacobian penalty
double MiNTEnergy::ComputeScaledJacobianPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeScaledJacobianPenaltyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the scaled jacobian penalty per tet
// I.e. 1 - det(f) / \prod |f_i| =
//      1 - det(f)/ \prod |f_i|
// IT IS PREFERABLE TO ENFORCE F TO HAVE POSITIVE DET EVERYWHERE
// double MiNTEnergy::ComputeScaledJacobianPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd*
// deriv,
//                                                       Eigen::MatrixXd* hess, bool is_PSD_proj) {
//     assert(frames.cols() == 9);

//     double energy = 0;
//     double weight = tet_volumes_[tet_id] * (1.0 / (ave_edge_len_ * ave_edge_len_)) *
//     this->weights_.w_scaled_jacobian;

//     Eigen::VectorXd deriv_det;
//     Eigen::Matrix<double, 9, 9> hess_det;
//     deriv_det.resize(9);
//     hess_det.setZero(9, 9);

//     if (deriv) {
//         deriv->resize(9);
//     }
//     if (hess) {
//         hess->setZero(9, 9);
//     }

//     Eigen::Matrix3d f;
//     f.row(0) = frames.row(tet_id).segment<3>(0);
//     f.row(1) = frames.row(tet_id).segment<3>(3);
//     f.row(2) = frames.row(tet_id).segment<3>(6);

//     double det = Det3x3(f, (deriv || hess) ? &deriv_det : nullptr, hess ? &hess_det : nullptr);

//     double f0n = f.row(0).norm();
//     double f1n = f.row(1).norm();
//     double f2n = f.row(2).norm();

//     double prod =  f0n * f1n * f2n + 1e-12;
//     double prod_inv = 1.0 / prod;

//     // energy = (1 - det * prod_inv) * weight;
//     // energy = prod_inv * weight;
//     energy = prod * weight;

//     // Compute derivatives if requested
//     Eigen::VectorXd d_prod(9), d_prod_inverse(9);
//     if (deriv) {
//         d_prod.setZero(9);
//         d_prod.segment<3>(0) = (f1n * f2n / f0n) * f.row(0);
//         d_prod.segment<3>(3) = (f0n * f2n / f1n) * f.row(1);
//         d_prod.segment<3>(6) = (f0n * f1n / f2n) * f.row(2);

//         d_prod_inverse = -prod_inv * prod_inv * d_prod;

//         // *deriv = weight * (-det * d_prod_inverse + prod_inv * deriv_det);
//         *deriv = weight * ( d_prod );

//     }

//     // Compute Hessians if requested
//     if (hess) {
//         Eigen::Matrix<double, 9, 9> d_d_prod, d_d_prod_inverse;

//         d_d_prod.setZero(9, 9);

//         Eigen::Matrix<double, 9, 9> H; // Hessian matrix

//         // Diagonal blocks
//         d_d_prod.block<3, 3>(0, 0) = (f1n * f2n / (f0n * f0n* f0n)) * Eigen::Matrix3d::Identity();
//         d_d_prod.block<3, 3>(3, 3) = (f0n * f2n / (f1n * f1n* f1n)) * Eigen::Matrix3d::Identity();
//         d_d_prod.block<3, 3>(6, 6) = (f0n * f1n / (f2n * f2n* f2n)) * Eigen::Matrix3d::Identity();

//         // Off-diagonal blocks
//         d_d_prod.block<3, 3>(0, 3) = -(f2n / (f0n * f1n)) * (f.row(0).transpose() * f.row(1));
//         d_d_prod.block<3, 3>(0, 6) = -(f1n / (f0n * f2n)) * (f.row(0).transpose() * f.row(2));
//         d_d_prod.block<3, 3>(3, 0) = d_d_prod.block<3, 3>(0, 3).transpose();
//         d_d_prod.block<3, 3>(3, 6) = -(f0n / (f1n * f2n)) * (f.row(1).transpose() * f.row(2));
//         d_d_prod.block<3, 3>(6, 0) = d_d_prod.block<3, 3>(0, 6).transpose();
//         d_d_prod.block<3, 3>(6, 3) = d_d_prod.block<3, 3>(3, 6).transpose();

//         // Compute d_d_prod_inverse
//         d_d_prod_inverse = -prod_inv * prod_inv * d_d_prod + 2.0 * prod_inv * prod_inv * prod_inv * d_prod *
//         d_prod.transpose();

//         // *hess = weight * (d_d_prod_inverse * (-det) + d_prod_inverse * deriv_det.transpose() + deriv_det *
//         d_prod_inverse.transpose() + prod_inv * hess_det);

//         *hess = weight * (H );

//         if (is_PSD_proj) {
//             *hess = SelfAdjSPDProjection(*hess);
//         }
//     }

//     return energy;
// }

// Compute the scaled jacobian penalty per tet
// I.e. 1 - pow(det(f) / \prod |f_i|, 2) =
//      1 - pow(det(f), 2)/ \prod |f_i|^2
double MiNTEnergy::ComputeScaledJacobianPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                      Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];

    weight *= (1.0 / ave_edge_len_ / ave_edge_len_);   // make it unit of length, scaled by coeff

    // if (tet_id > mesh_not_extended_->nTets())
    // {
    //     weight *= this->weights_.w_bound / 5.;
    //     weight *= this->weights_.w_bound_rescale;

    // }
    // else {
    // weight *= this->weights_.w_scaled_jacobian;

    // }

    // if (tet_id > mesh_not_extended_->nTets())
    // {
    //     weight *= this->weights_.w_scaled_jacobian * this->weights_.w_bound_rescale;
    // }
    // else {

    //     // weight *= this->weights_.w_smooth * 5;
    //     weight *= this->weights_.w_scaled_jacobian;

    // }

    if (tet_id > mesh_not_extended_->nTets()) {
        weight *= this->weights_.w_smooth * this->weights_.w_bound_rescale * 1e-2;   //  .001; //
    } else {
        // weight *= this->weights_.w_smooth * 5;
        weight *= this->weights_.w_smooth * 5;   // * .001; //
    }
    // weight = 0;

    Eigen::VectorXd deriv_det;
    Eigen::Matrix<double, 9, 9> hess_det;
    deriv_det.resize(9);
    hess_det.setZero(9, 9);

    if (deriv) {
        deriv->resize(9);
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    double det = Det3x3(f, (deriv || hess) ? &deriv_det : nullptr, hess ? &hess_det : nullptr);

    double prod_sq = f.row(0).squaredNorm() * f.row(1).squaredNorm() * f.row(2).squaredNorm() + 1e-12;
    double det_sq = det * det;
    double prod_sq_inv = 1. / prod_sq;

    energy = 1 - det_sq * prod_sq_inv;
    energy *= weight;

    // compute primatives

    Eigen::VectorXd d_det_squared;
    Eigen::MatrixXd d_d_det_squared;

    Eigen::VectorXd d_prod_sq;
    Eigen::MatrixXd d_d_prod_sq;

    Eigen::VectorXd d_prod_sq_inverse;
    Eigen::MatrixXd d_d_prod_sq_inverse;
    if (deriv) {
        // d: energy = 1 - pow(det, 2);
        d_det_squared = -2. * det * (deriv_det);
        d_prod_sq.setZero(9);
        d_prod_sq.segment<3>(0) = 2 * f.row(0) * f.row(1).squaredNorm() * f.row(2).squaredNorm();
        d_prod_sq.segment<3>(3) = 2 * f.row(1) * f.row(0).squaredNorm() * f.row(2).squaredNorm();
        d_prod_sq.segment<3>(6) = 2 * f.row(2) * f.row(0).squaredNorm() * f.row(1).squaredNorm();

        d_prod_sq_inverse = -1. / (prod_sq * prod_sq) * d_prod_sq;
    }
    if (hess) {
        d_d_det_squared = -2 * deriv_det * deriv_det.transpose() - 2 * det * hess_det;

        d_d_prod_sq.setZero(9, 9);
        d_d_prod_sq.block<3, 3>(0, 0) =
            2 * f.row(1).squaredNorm() * f.row(2).squaredNorm() * Eigen::Matrix3d::Identity();
        d_d_prod_sq.block<3, 3>(3, 3) =
            2 * f.row(0).squaredNorm() * f.row(2).squaredNorm() * Eigen::Matrix3d::Identity();
        d_d_prod_sq.block<3, 3>(6, 6) =
            2 * f.row(0).squaredNorm() * f.row(1).squaredNorm() * Eigen::Matrix3d::Identity();

        d_d_prod_sq.block<3, 3>(0, 6) = 4 * f.row(0).transpose() * f.row(2) * f.row(1).squaredNorm();
        d_d_prod_sq.block<3, 3>(0, 3) = 4 * f.row(0).transpose() * f.row(1) * f.row(2).squaredNorm();
        d_d_prod_sq.block<3, 3>(3, 0) = 4 * f.row(1).transpose() * f.row(0) * f.row(2).squaredNorm();
        d_d_prod_sq.block<3, 3>(3, 6) = 4 * f.row(1).transpose() * f.row(2) * f.row(0).squaredNorm();
        d_d_prod_sq.block<3, 3>(6, 0) = 4 * f.row(2).transpose() * f.row(0) * f.row(1).squaredNorm();
        d_d_prod_sq.block<3, 3>(6, 3) = 4 * f.row(2).transpose() * f.row(1) * f.row(0).squaredNorm();

        d_d_prod_sq_inverse = -1. / (prod_sq * prod_sq) * d_d_prod_sq +
                              2. / (prod_sq * prod_sq * prod_sq) * d_prod_sq * d_prod_sq.transpose();
    }

    if (deriv) {
        (*deriv) = weight * (-det_sq * d_prod_sq_inverse + (1. / prod_sq) * d_det_squared);
    }
    if (hess) {
        (*hess) = -det_sq * d_d_prod_sq_inverse + d_det_squared * d_prod_sq_inverse.transpose();
        (*hess) += prod_sq_inv * d_d_det_squared + d_prod_sq_inverse * d_det_squared.transpose();
        (*hess) *= weight;
    }

    if (hess && is_PSD_proj) {
        (*hess) = SelfAdjSPDProjection(*hess);
    }
    return energy;
}

// stable neohookean for smaller poissons ratio with barrier from degenerating to a point
double MiNTEnergy::ComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                   Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];
    double base_weight = weight * (1.0 / ave_edge_len_ / ave_edge_len_);
    // weight = base_weight * this->weights_.w_smooth * 100000.;
    weight = base_weight * this->weights_.w_smooth * 1.;
    if (tet_id > mesh_not_extended_->nTets()) {
        weight = 0.;
        // weight = base_weight * this->weights_.w_smooth * .001 * this->weights_.w_bound_rescale;
    }

    // weight = 0;

    Eigen::VectorXd deriv_det(9);
    Eigen::Matrix<double, 9, 9> hess_det;
    deriv_det.setZero();
    hess_det.setZero();

    if (deriv) {
        deriv->resize(9);
        deriv->setZero();
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    double det = Det3x3(f, (deriv || hess) ? &deriv_det : nullptr, hess ? &hess_det : nullptr);
    // double mu = 1.;
    // double mu = 1e-2;

    // double mu = 0.2;
    // double mu = 8. / 5. - 1e-4;
    double mu = 8. / 5. - 1e-1;

    // double mu = 7. / 5. - 1e-4;
    // double mu = 5. / 5.;
    // (1 − 5/8 * (7/5)) / (2*(1 + 1/8  (7/5))) <--- poisson ratio

    double lambda = 1.0;
    double alpha = 1 + mu / lambda - mu / (4 * lambda);
    double I_C = f.row(0).squaredNorm() + f.row(1).squaredNorm() + f.row(2).squaredNorm();
    double det_diff = det - alpha;

    // Energy calculation
    energy += mu / 2 * (I_C - 3) * weight;
    energy += lambda / 2 * det_diff * det_diff * weight;
    energy -= mu / 2 * log(I_C + 1) * weight;

    if (deriv) {
        Eigen::VectorXd grad_energy = Eigen::VectorXd::Zero(9);

        // Gradient of the determinant term
        grad_energy += lambda * det_diff * deriv_det * weight;

        // Gradient of the squared norm term
        for (int i = 0; i < 3; ++i) {
            Eigen::Vector3d v = f.row(i);
            Eigen::Vector3d grad_diff = mu * weight * v;
            grad_energy.segment<3>(3 * i) += grad_diff;
        }

        // Gradient of the logarithm term
        for (int i = 0; i < 3; ++i) {
            Eigen::Vector3d v = f.row(i);
            Eigen::Vector3d grad_log = (mu * weight / (I_C + 1)) * v;
            grad_energy.segment<3>(3 * i) -= grad_log;
        }

        *deriv = grad_energy;
    }

    if (hess) {
        Eigen::Matrix<double, 9, 9> hess_energy = Eigen::Matrix<double, 9, 9>::Zero();

        // Hessian of the determinant term
        hess_energy += lambda * weight * (deriv_det * deriv_det.transpose() + det_diff * hess_det);

        // Hessian of the squared norm term
        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix3d hess_diff = mu * weight * Eigen::Matrix3d::Identity();
            hess_energy.block<3, 3>(3 * i, 3 * i) += hess_diff;
        }

        // Hessian of the logarithm term
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i == j) {
                    Eigen::Vector3d v = f.row(i);
                    Eigen::Matrix3d hess_log = (mu * weight / (I_C + 1)) * Eigen::Matrix3d::Identity() -
                                               (2.0 * mu * weight / pow(I_C + 1, 2)) * (v * v.transpose());
                    hess_energy.block<3, 3>(3 * i, 3 * i) -= hess_log;
                } else {
                    Eigen::Vector3d u = f.row(i);
                    Eigen::Vector3d v = f.row(j);

                    Eigen::Matrix3d hess_log = -(mu * weight / pow(I_C + 1, 2)) * (u * v.transpose());
                    hess_energy.block<3, 3>(3 * i, 3 * j) -= hess_log;
                    hess_energy.block<3, 3>(3 * j, 3 * i) -= hess_log.transpose();   // Symmetric contribution
                }
            }
        }

        // hess_energy += (mu * weight / (I_C + 1)) * Eigen::Matrix<double, 9, 9>::Identity();

        // for (int i = 0; i < 3; ++i) {
        //     for (int j = 0; j < 3; ++j) {
        //         int k = 3 * i + j; // Index for the gradient element
        //         hess_energy(k, k) -= 2 * (mu * weight / std::pow(I_C + 1, 2)) * f(i, j) * f(i, j);
        //     }
        // }

        // // // Hessian of the logarithm term (Corrected)
        // for (int i = 0; i < 3; ++i) {
        //     Eigen::Vector3d v = f.row(i);
        //     Eigen::Matrix3d hess_log = -2 * (mu * weight / pow(I_C + 1, 2)) * v * v.transpose();  // Only keep the
        //     second derivative part hess_energy.block<3, 3>(3 * i, 3 * i) += hess_log;
        // }

        *hess = hess_energy;
    }

    if (hess && is_PSD_proj) {
        *hess = SelfAdjSPDProjection(*hess);
    }

    return energy;
}

// ComputeARAPEnergyPerTet

// double MiNTEnergy::ComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
//                                            Eigen::MatrixXd* hess, bool is_PSD_proj) {
//     assert(frames.cols() == 9);
//     double energy = 0;
//     double weight = tet_volumes_[tet_id];
//     double base_weight = weight * (1.0 / ave_edge_len_ / ave_edge_len_);

//     double eps_weight = 0; // hack...
//     weight = base_weight * this->weights_.w_smooth * .001;
//     if (tet_id > mesh_not_extended_->nTets()) {
//         eps_weight = 0;
//         weight = 0.;
//     }

//     Eigen::VectorXd deriv_arap(9);
//     Eigen::Matrix<double, 9, 9> hess_arap;
//     deriv_arap.setZero();
//     hess_arap.setZero();

//     if (deriv) {
//         deriv->resize(9);
//         deriv->setZero();
//     }
//     if (hess) {
//         hess->setZero(9, 9);
//     }

//     Eigen::Matrix3d f;
//     f.row(0) = frames.row(tet_id).segment<3>(0);
//     f.row(1) = frames.row(tet_id).segment<3>(3);
//     f.row(2) = frames.row(tet_id).segment<3>(6);

//     Eigen::Matrix3d cauchy_green = f.transpose() * f;
//     Eigen::JacobiSVD<Eigen::MatrixXd> svd(cauchy_green, Eigen::ComputeFullU | Eigen::ComputeFullV);
//     Eigen::Matrix3d r = svd.matrixU() * svd.matrixV().transpose();

//     Eigen::Matrix3d diff = f - r;
//     double frob_norm = diff.squaredNorm();

//     // Energy calculation
//     energy = 0.5 * weight * frob_norm;

//     if (deriv) {
//         Eigen::Matrix3d grad_energy = weight * diff;

//         for (int i = 0; i < 3; ++i) {
//             deriv->segment<3>(3 * i) = grad_energy.row(i);
//         }
//     }

//     if (hess) {
//         Eigen::Matrix<double, 9, 9> I = Eigen::Matrix<double, 9, 9>::Identity();
//         Eigen::Matrix<double, 9, 9> hess_diff = I * weight;

//         for (int i = 0; i < 3; ++i) {
//             hess_diff.block<3, 3>(3 * i, 3 * i) -= weight * r.row(i).transpose() * r.row(i);
//         }

//         *hess = hess_diff;
//     }

//     if (hess && is_PSD_proj) {
//         *hess = SelfAdjSPDProjection(*hess);
//     }

//     return energy;
// }

// ComputeARAPEnergyPerTet

double MiNTEnergy::ComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];
    double base_weight = weight * (1.0 / ave_edge_len_ / ave_edge_len_);
    double mu = this->weights_.w_smooth; // assuming this->weights_.w_smooth is the material constant mu

    weight = base_weight * 0.05;
    if (tet_id > mesh_not_extended_->nTets()) {
        weight = 0.;
    }

    if (deriv) {
        deriv->resize(9);
        deriv->setZero();
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d F;
    F.row(0) = frames.row(tet_id).segment<3>(0);
    F.row(1) = frames.row(tet_id).segment<3>(3);
    F.row(2) = frames.row(tet_id).segment<3>(6);

    // Compute R using polar decomposition F = RS
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();

    // Compute energy
    double normF = F.squaredNorm();
    double trFR = (F.transpose() * R).trace();
    energy = (mu / 2) * (normF + 3 - 2 * trFR) * weight;

    if (deriv) {
        // Compute gradient
        Eigen::Matrix3d grad_matrix = mu * (F - R) * weight;
        deriv->segment<3>(0) = grad_matrix.row(0);
        deriv->segment<3>(3) = grad_matrix.row(1);
        deriv->segment<3>(6) = grad_matrix.row(2);
    }

    if (hess) {
        // Compute Hessian
        Eigen::Matrix<double, 9, 9> hess_energy;
        hess_energy.setZero();
        Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();

        for (int i = 0; i < 3; ++i) {
            hess_energy.block<3, 3>(3 * i, 3 * i) = mu * I3 * weight;
        }

        *hess = hess_energy;
    }

    if (hess && is_PSD_proj) {
        *hess = SelfAdjSPDProjection(*hess);
    }

    return energy;
}



// bonnet and wood  ComputeBonetWoodEnergy
// Function to compute Bonet and Wood (2008) version of Neo-Hookean energy per tetrahedron


// stable neohookean
double MiNTEnergy::ComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                   Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];
    double base_weight = weight * (1.0 / ave_edge_len_ / ave_edge_len_);

    double eps_weight = 0; // hack...
    weight = base_weight * this->weights_.w_smooth * .1;
    if (tet_id > mesh_not_extended_->nTets()) {
        // eps_weight  = base_weight * this->weights_.w_smooth * 1e-8;

        eps_weight =  0; // base_weight * this->weights_.w_smooth * .0001;
        // weight  = base_weight * this->weights_.w_smooth * .0001;

        weight = 0.;
    } else {

    }
    // weight = 0;
    // weight = 1;

    Eigen::VectorXd deriv_det(9);
    Eigen::Matrix<double, 9, 9> hess_det;
    deriv_det.setZero();
    hess_det.setZero();

    if (deriv) {
        deriv->resize(9);
        deriv->setZero();
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    double det = Det3x3(f, (deriv || hess) ? &deriv_det : nullptr, hess ? &hess_det : nullptr);

    // double mu = .2; // ratio .4
    // double mu = 8./5. - 1e-6; // ratio 0
    // double mu = 4./5.; // ratio .22
    double mu = 1e-6;
    double lambda = 1.0;

    // if (tet_id > mesh_not_extended_->nTets()) {
    //     mu = 8./5. - 1e-6;
    //     lambda = 1.0;
    // }
    // double alpha = 1 + mu / lambda - mu / (4 * lambda); // this is for the log barrier term.
    double alpha = 1;
    double det_diff = det - alpha;

    // Energy calculation
    energy = (-mu * det_diff + 0.5 * lambda * det_diff * det_diff) * weight;
    energy += mu / 2 * (f.row(0).squaredNorm() + f.row(1).squaredNorm() + f.row(2).squaredNorm() - 3) * (weight +
eps_weight);


    if (deriv) {
        Eigen::VectorXd grad_energy = (-mu + lambda * det_diff) * deriv_det * weight;

        for (int i = 0; i < 3; ++i) {
            Eigen::Vector3d v = f.row(i);
            Eigen::Vector3d grad_diff = mu * (weight + eps_weight) * v;   // Gradient of squared norm term
            grad_energy.segment<3>(3 * i) += grad_diff;
        }

        *deriv = grad_energy;
    }

    if (hess) {
        Eigen::Matrix<double, 9, 9> hess_energy =
            lambda * weight * (deriv_det * deriv_det.transpose() + det_diff * hess_det) - mu * weight * hess_det;

        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix3d hess_diff = mu * (weight + eps_weight) * Eigen::Matrix3d::Identity();   // Hessian of
squared norm term hess_energy.block<3, 3>(3 * i, 3 * i) += hess_diff;
        }

        *hess = hess_energy;
    }

    if (hess && is_PSD_proj) {
        *hess = SelfAdjSPDProjection(*hess);
    }

    return energy;
}



///////////////////////////////////

// Compute the scaled jacobian penalty per tet
// I.e. \prod |f_i|^2 / pow(det(f), 2) - 1
// like the scaled jacobian, this energy is also zero when frames are orthogonal,
// but it becomes unbounded when frames degenerate.
// maybe useful to add in a "symmetric" term like this?
double MiNTEnergy::ComputeScaledJacobianInversePenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                                             Eigen::VectorXd* deriv, Eigen::MatrixXd* hess,
                                                             bool is_PSD_proj) {
    assert(frames.cols() == 9);

    double energy = 0;
    double weight = tet_volumes_[tet_id];
    weight *= (1.0 / (ave_edge_len_ * ave_edge_len_));   // make it unit of length, scaled by coeff
    weight *= this->weights_.w_scaled_jacobian;

    Eigen::VectorXd deriv_det(9);
    Eigen::Matrix<double, 9, 9> hess_det;
    hess_det.setZero();

    if (deriv) {
        deriv->resize(9);
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    double det = Det3x3(f, (deriv || hess) ? &deriv_det : nullptr, hess ? &hess_det : nullptr);

    double prod_sq = f.row(0).squaredNorm() * f.row(1).squaredNorm() * f.row(2).squaredNorm() + 1e-12;
    double det_sq = det * det;
    double det_sq_inv = 1. / (det_sq + 1e-12);

    energy = prod_sq * det_sq_inv - 1;
    energy *= weight;

    // Compute derivatives
    Eigen::VectorXd d_det_squared;
    Eigen::VectorXd d_det_sq_inverse;
    Eigen::VectorXd d_prod_sq(9);
    if (deriv) {
        d_det_squared = 2. * det * deriv_det;
        d_det_sq_inverse = -2.0 * det * det_sq_inv * det_sq_inv * deriv_det;

        d_prod_sq.setZero();
        d_prod_sq.segment<3>(0) = 2 * f.row(0) * f.row(1).squaredNorm() * f.row(2).squaredNorm();
        d_prod_sq.segment<3>(3) = 2 * f.row(1) * f.row(0).squaredNorm() * f.row(2).squaredNorm();
        d_prod_sq.segment<3>(6) = 2 * f.row(2) * f.row(0).squaredNorm() * f.row(1).squaredNorm();
    }

    // Compute Hessians
    Eigen::MatrixXd d_d_det_squared;
    Eigen::MatrixXd d_d_det_sq_inverse;
    Eigen::MatrixXd d_d_prod_sq(9, 9);

    if (hess) {
        d_d_det_squared = 2 * deriv_det * deriv_det.transpose() + 2 * det * hess_det;
        d_d_det_sq_inverse = -det_sq_inv * det_sq_inv * d_d_det_squared +
                             2. * det_sq_inv * det_sq_inv * det_sq_inv * d_det_squared * d_det_squared.transpose();

        d_d_prod_sq.setZero();
        d_d_prod_sq.block<3, 3>(0, 0) =
            2 * f.row(1).squaredNorm() * f.row(2).squaredNorm() * Eigen::Matrix3d::Identity();
        d_d_prod_sq.block<3, 3>(3, 3) =
            2 * f.row(0).squaredNorm() * f.row(2).squaredNorm() * Eigen::Matrix3d::Identity();
        d_d_prod_sq.block<3, 3>(6, 6) =
            2 * f.row(0).squaredNorm() * f.row(1).squaredNorm() * Eigen::Matrix3d::Identity();

        d_d_prod_sq.block<3, 3>(0, 6) = 4 * f.row(0).transpose() * f.row(2) * f.row(1).squaredNorm();
        d_d_prod_sq.block<3, 3>(0, 3) = 4 * f.row(0).transpose() * f.row(1) * f.row(2).squaredNorm();
        d_d_prod_sq.block<3, 3>(3, 0) = 4 * f.row(1).transpose() * f.row(0) * f.row(2).squaredNorm();
        d_d_prod_sq.block<3, 3>(3, 6) = 4 * f.row(1).transpose() * f.row(2) * f.row(0).squaredNorm();
        d_d_prod_sq.block<3, 3>(6, 0) = 4 * f.row(2).transpose() * f.row(0) * f.row(1).squaredNorm();
        d_d_prod_sq.block<3, 3>(6, 3) = 4 * f.row(2).transpose() * f.row(1) * f.row(0).squaredNorm();
    }

    if (deriv) {
        *deriv = weight * (prod_sq * d_det_sq_inverse + det_sq_inv * d_prod_sq);
    }
    if (hess) {
        *hess = prod_sq * d_d_det_sq_inverse + d_prod_sq * d_det_sq_inverse.transpose() + det_sq_inv * d_d_prod_sq +
                d_det_sq_inverse * d_prod_sq.transpose();
        *hess *= weight;
    }

    if (hess && is_PSD_proj) {
        *hess = SelfAdjSPDProjection(*hess);
    }

    return energy;
}

///////////// Unit Barrier Energy //////////////
// Compute the unit barrier penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
// frames per tetrahedron
double MiNTEnergy::ComputeUnitBarrierPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                                         std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                         bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int tet_id, Eigen::VectorXd* tet_deriv,
                           Eigen::MatrixXd* tet_hess, bool is_tet_PSD_proj) {
        return ComputeUnitBarrierPenaltyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
        // return ComputeCrossProductEnergyPerTet(input_frames, tet_id, tet_deriv, tet_hess, is_tet_PSD_proj);
    };
    return AssembleEnergyFromTets(frames, deriv, hess_triplets, is_PSD_proj, energy_func);
}

// Compute the scaled jacobian penalty
double MiNTEnergy::ComputeUnitBarrierPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv,
                                             Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
    std::vector<Eigen::Triplet<double>> T;
    double energy = ComputeUnitBarrierPenaltyWithTriplets(frames, deriv, hess ? &T : nullptr, is_PSD_proj);

    if (hess) {
        hess->resize(frames.rows() * frames.cols(), frames.rows() * frames.cols());
        hess->setFromTriplets(T.begin(), T.end());
    }

    return energy;
}

// Compute the scaled jacobian penalty per tet
// I.e. 1 - pow(det(f) / \prod |f_i|, 2) =
//      1 - pow(det(f), 2)/ \prod |f_i|^2
double MiNTEnergy::ComputeUnitBarrierPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                   Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];
    weight *= (this->weights_.w_unit_barrier / ave_edge_len_ / ave_edge_len_);

    Eigen::VectorXd deriv_det;
    Eigen::Matrix<double, 9, 9> hess_det;
    deriv_det.resize(9);
    hess_det.setZero(9, 9);

    if (deriv) {
        deriv->setZero(9);
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    double max_norm = f.row(0).norm();
    int max_idx = 0;

    double tie_break_thresh = 1e-6;

    if (f.row(1).norm() > max_norm + tie_break_thresh) {
        max_norm = f.row(1).norm();
        max_idx = 1;
    } else if (f.row(2).norm() > max_norm + tie_break_thresh) {
        max_norm = f.row(2).norm();
        max_idx = 2;
    }

    max_norm = max_norm + 1e-12;

    energy = pow(1 - 1 / (max_norm), 2);
    energy *= weight;

    // std::cout << "max_norm: " << max_norm << " max_idx: " << max_idx << " energy: " << energy << std::endl;

    // compute primatives

    if (deriv) {
        // deriv->segment<3>(3 * max_idx) = -2. / (pow(max_norm, 3)) * f.row(max_idx) * (1 - 1 / (max_norm));
        // Calculate t and t1
        // double t = 1.0 / max_norm;
        // double t1 = 1.0 / (max_norm * max_norm * max_norm);

        // // Compute the scalar factor for the gradient
        // double scalar_factor = 2.0 * t1 * (1.0 - t);

        // Update the derivative vector using segment at 3 * max_idx
        // deriv->segment<3>(3 * max_idx) = scalar_factor * f.row(max_idx).transpose();

        //

        Eigen::Vector3d x = f.row(max_idx).transpose();

        deriv->segment<3>(3 * max_idx) =
            2.0 * 1. / std::pow((x.transpose() * x), 1.5) * (1.0 - 1. / std::pow(x.transpose() * x, .5)) * x;
        (*deriv) *= weight;
        // deriv->segment<3>(3 * max_idx) = 2*( std::pow(x.transpose()*x, (1/2)) ) *
        // (1-std::pow((x.transpose()*x),(-1/2)))*x;
    }

    if (hess) {
        // via matrixcalculus.org
        // input: 2.0 * 1/((x'*x)^(3/2)) * (1.0 - 1/((x'*x)^(1/2))) * x

        Eigen::Vector3d x = f.row(max_idx).transpose();

        // Compute t0 = x^T x (norm squared of x)
        double t0 = x.squaredNorm();

        // Compute t and other coefficients
        double t = std::sqrt(t0);              // t = \|x\|
        double t1 = 1.5;                       // t1 = 3/2
        double t2 = 0.5;                       // t2 = 1/2
        double t3 = 1.0 - std::pow(t0, -t2);   // t3 = 1.0 - (x^T x)^{-t2}
        double t4 = std::pow(t0, -t1);         // t4 = (x^T x)^{-t1}

        // Compute outer product T5 = x x^T
        Eigen::MatrixXd T5 = x * x.transpose();

        // Compute identity matrix I
        Eigen::MatrixXd I = Eigen::Matrix3d::Identity();

        // Compute Hessian terms
        Eigen::MatrixXd HessianTerm1 = 2.0 * t3 * t4 * I;
        Eigen::MatrixXd HessianTerm2 = 6.0 * t3 * std::pow(t0, (-(1.0 + t1))) * T5;
        Eigen::MatrixXd HessianTerm3 = 2.0 * t4 * std::pow(t0, (-(1.0 + t2))) * T5;

        // Compute Hessian matrix H = HessianTerm1 - (HessianTerm2 + HessianTerm3)
        hess->block<3, 3>(3 * max_idx, 3 * max_idx) = HessianTerm1 - HessianTerm2 + HessianTerm3;

        (*hess) *= weight;
    }

    if (hess && is_PSD_proj) {
        hess->block<3, 3>(3 * max_idx, 3 * max_idx) = SelfAdjSPDProjection(hess->block<3, 3>(3 * max_idx, 3 * max_idx));
        // }
    }
    return energy;
}

double MiNTEnergy::ComputeCrossProductEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv,
                                                   Eigen::MatrixXd* hess, bool is_PSD_proj) {
    assert(frames.cols() == 9);
    double energy = 0;
    double weight = tet_volumes_[tet_id];
    weight *= (this->weights_.w_unit_barrier / ave_edge_len_ / ave_edge_len_);

    if (deriv) {
        deriv->setZero(9);
    }
    if (hess) {
        hess->setZero(9, 9);
    }

    Eigen::Matrix3d f;
    f.row(0) = frames.row(tet_id).segment<3>(0);
    f.row(1) = frames.row(tet_id).segment<3>(3);
    f.row(2) = frames.row(tet_id).segment<3>(6);

    auto compute_term = [&](const Eigen::Vector3d& a, const Eigen::Vector3d& b, int a_idx, int b_idx,
                            Eigen::VectorXd* deriv, Eigen::MatrixXd* hess) -> double {
        Eigen::Vector3d cross_product = a.cross(b);
        double norm_sq = cross_product.squaredNorm();
        double term = 1 - norm_sq;
        double energy_term = term * term;

        if (deriv) {
            Eigen::Matrix3d da_cross_db;
            da_cross_db << 0, -b(2), b(1), b(2), 0, -b(0), -b(1), b(0), 0;

            Eigen::Matrix3d db_cross_da;
            db_cross_da << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;

            Eigen::Vector3d grad_a = -4 * term * da_cross_db * cross_product;
            Eigen::Vector3d grad_b = 4 * term * db_cross_da * cross_product;

            deriv->segment<3>(a_idx) += grad_a;
            deriv->segment<3>(b_idx) += grad_b;

            if (hess) {
                Eigen::Matrix3d da_cross_db_outer =
                    da_cross_db * cross_product * cross_product.transpose() * da_cross_db.transpose();
                Eigen::Matrix3d db_cross_da_outer =
                    db_cross_da * cross_product * cross_product.transpose() * db_cross_da.transpose();
                Eigen::Matrix3d mixed_term =
                    da_cross_db * cross_product * cross_product.transpose() * db_cross_da.transpose();

                hess->block<3, 3>(a_idx, a_idx) +=
                    8 * (cross_product.squaredNorm() * da_cross_db * da_cross_db.transpose() - da_cross_db_outer);
                hess->block<3, 3>(b_idx, b_idx) +=
                    8 * (cross_product.squaredNorm() * db_cross_da * db_cross_da.transpose() - db_cross_da_outer);
                hess->block<3, 3>(a_idx, b_idx) -=
                    8 * (cross_product.squaredNorm() * da_cross_db * db_cross_da.transpose() - mixed_term);
                hess->block<3, 3>(b_idx, a_idx) = hess->block<3, 3>(a_idx, b_idx).transpose();
            }
        }

        return energy_term;
    };

    energy += compute_term(f.row(0), f.row(1), 0, 3, deriv, hess);
    energy += compute_term(f.row(0), f.row(2), 0, 6, deriv, hess);
    energy += compute_term(f.row(1), f.row(2), 3, 6, deriv, hess);

    energy *= weight;

    if (deriv) {
        (*deriv) *= weight;
    }

    if (hess) {
        (*hess) *= weight;
        if (is_PSD_proj) {
            (*hess) = SelfAdjSPDProjection(*hess);
        }
    }

    return energy;
}

///////////// Test functions //////////////
// Test function for general function
void MiNTEnergy::TestGeneralFunction(
    const Eigen::MatrixXd& frames,
    std::function<double(const Eigen::MatrixXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> energy_func) {
    Eigen::VectorXd deriv;
    Eigen::SparseMatrix<double> hess;

    double energy = energy_func(frames, &deriv, &hess, false);

    Eigen::MatrixXd frames1 = frames;

    Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(deriv.size());

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(10, -i);

        for (int j = 0; j < mesh_->nTets(); j++) {
            frames1.row(j) = frames.row(j) + eps * perturb_vec.segment<9>(9 * j).transpose();
        }

        Eigen::VectorXd deriv1;
        double energy1 = energy_func(frames1, &deriv1, nullptr, false);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
        std::cout << "deriv-hess: " << ((deriv1 - deriv) / eps - hess * perturb_vec).norm() << std::endl;
    }
}

// Test function for general function defined per face
void MiNTEnergy::TestGeneralFunctionPerFace(
    const Eigen::MatrixXd& frames, int face_id,
    std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func) {
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);

    if (tet_id0 == -1 || tet_id1 == -1) {
        std::cout << "Boundary face" << std::endl;
        return;
    }

    Eigen::VectorXd deriv;
    Eigen::MatrixXd hess;

    double energy = energy_func(frames, face_id, &deriv, &hess, false);

    Eigen::MatrixXd frames1 = frames;

    Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(18);

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(10, -i);
        frames1.row(tet_id0) = frames.row(tet_id0) + eps * perturb_vec.segment<9>(0).transpose();
        frames1.row(tet_id1) = frames.row(tet_id1) + eps * perturb_vec.segment<9>(9).transpose();

        Eigen::VectorXd deriv1;
        double energy1 = energy_func(frames1, face_id, &deriv1, nullptr, false);

        std::cout << "eps: " << eps << std::endl;
        std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
        std::cout << "deriv-hess: " << ((deriv1 - deriv) / eps - hess * perturb_vec).norm() << std::endl;
    }
}

// Test function for general function defined per tet
void MiNTEnergy::TestGeneralFunctionPerTet(
    const Eigen::MatrixXd& frames, int tet_id,
    std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func) {
    Eigen::VectorXd deriv;
    Eigen::MatrixXd hess;

    double energy = energy_func(frames, tet_id, &deriv, &hess, false);

    Eigen::MatrixXd frames1 = frames;

    Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(9);

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(10, -i);
        frames1.row(tet_id) = frames.row(tet_id) + eps * perturb_vec.transpose();

        Eigen::VectorXd deriv1;
        double energy1 = energy_func(frames1, tet_id, &deriv1, nullptr, false);

        double hess_err = ((deriv1 - deriv) / eps - hess * perturb_vec).norm();

        std::cout << "eps: " << eps << std::endl;
        std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
        std::cout << "deriv-hess: " << hess_err << std::endl;

        if (i == 6 && hess_err > 1e-2) {
            Eigen::Matrix<double, 9, 9> hess_debug;
            Eigen::Matrix<double, 9, 9> hess_finite_diff;
            for (int j = 0; j < 9; j++) {
                Eigen::VectorXd pvec = Eigen::VectorXd::Zero(9);
                pvec(j) = 1.;
                frames1.row(tet_id) = frames.row(tet_id) + eps * pvec.transpose();
                Eigen::VectorXd deriv1;
                double energy1 = energy_func(frames1, tet_id, &deriv1, nullptr, false);
                Eigen::VectorXd hess_col = ((deriv1 - deriv) / eps - hess * pvec);
                for (int k = 0; k < 9; k++) {
                    if (hess_col(k) < 1e-8) {
                        hess_col(k) = 0;
                    }
                }
                hess_debug.row(j) = hess_col;
                hess_finite_diff.row(j) = (deriv1 - deriv) / eps;
            }
            std::cout << "finite diff hessian" << std::endl << hess_finite_diff << std::endl;
            std::cout << "true hessian" << std::endl << hess << std::endl;
            std::cout << "round differences" << std::endl << hess_debug << std::endl;
        }
    }
}

// Test function for computing Moment tensors per tet
void MiNTEnergy::TestComputeMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd frames1 = frames;
    for (int k = 2; k <= 2; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        Eigen::VectorXd moment_tensor;
        Eigen::MatrixXd tensor_deriv;
        std::vector<Eigen::MatrixXd> tensor_hess;

        moment_tensor = ComputeMomentTensorsPerTet(frames, tet_id, k, &tensor_deriv, &tensor_hess);

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(9);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id) = frames.row(tet_id) + eps * perturb_vec.transpose();

            Eigen::VectorXd moment_tensor1;
            Eigen::MatrixXd tensor_deriv1;

            moment_tensor1 = ComputeMomentTensorsPerTet(frames1, tet_id, k, &tensor_deriv1, nullptr);

            std::cout << "\neps: " << eps << std::endl;
            std::cout << "grad-check: " << ((moment_tensor1 - moment_tensor) / eps - tensor_deriv * perturb_vec).norm()
                      << std::endl;
            std::cout << "hess-check: " << std::endl;
            for (int j = 0; j < moment_tensor.size(); j++) {
                std::cout << j << "-th entry: "
                          << ((tensor_deriv1.row(j) - tensor_deriv.row(j)) / eps -
                              (tensor_hess[j] * perturb_vec).transpose())
                                 .norm()
                          << std::endl;
            }
        }
    }
}

// Test function for computing Kruskal tensors per tet
void MiNTEnergy::TestComputeKruskalTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd frames1 = frames;
    for (int k = 2; k <= 2; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        Eigen::VectorXd kruskal_tensor;
        Eigen::MatrixXd tensor_deriv;
        std::vector<Eigen::MatrixXd> tensor_hess;

        kruskal_tensor = ComputeKruskalTensorsPerTet(frames, tet_id, k, &tensor_deriv, &tensor_hess);

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(9);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id) = frames.row(tet_id) + eps * perturb_vec.transpose();

            Eigen::VectorXd kruskal_tensor1;
            Eigen::MatrixXd tensor_deriv1;

            kruskal_tensor1 = ComputeKruskalTensorsPerTet(frames1, tet_id, k, &tensor_deriv1, nullptr);

            std::cout << "\neps: " << eps << std::endl;
            std::cout << "grad-check: "
                      << ((kruskal_tensor1 - kruskal_tensor) / eps - tensor_deriv * perturb_vec).norm() << std::endl;
            std::cout << "hess-check: " << std::endl;
            for (int j = 0; j < kruskal_tensor.size(); j++) {
                std::cout << j << "-th entry: "
                          << ((tensor_deriv1.row(j) - tensor_deriv.row(j)) / eps -
                              (tensor_hess[j] * perturb_vec).transpose())
                                 .norm()
                          << std::endl;
            }
        }
    }
}

// Test function for computing tensors per tet facet
void MiNTEnergy::TestComputeFaceProjectedTensors(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd frames1 = frames;
    for (int k = 2; k <= 2; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        Eigen::VectorXd kruskal_tensor;
        Eigen::MatrixXd tensor_deriv;
        std::vector<Eigen::MatrixXd> tensor_hess;

        bool is_kruskal = true;
        int facet_id = 1;

        kruskal_tensor = ComputeBasisProjectedTensorsPerTetFace(tet_facet_basis_.at(tet_id).at(facet_id), frames,
                                                                is_kruskal, tet_id, k, &tensor_deriv, &tensor_hess);

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(9);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id) = frames.row(tet_id) + eps * perturb_vec.transpose();

            Eigen::VectorXd kruskal_tensor1;
            Eigen::MatrixXd tensor_deriv1;

            // kruskal_tensor1 = ComputeKruskalTensorsPerTet(frames1, tet_id, k, &tensor_deriv1, nullptr);
            kruskal_tensor1 = ComputeBasisProjectedTensorsPerTetFace(tet_facet_basis_.at(tet_id).at(facet_id), frames1,
                                                                     is_kruskal, tet_id, k, &tensor_deriv1, nullptr);

            std::cout << "\neps: " << eps << std::endl;
            std::cout << "grad-check: "
                      << ((kruskal_tensor1 - kruskal_tensor) / eps - tensor_deriv * perturb_vec).norm() << std::endl;
            std::cout << "hess-check: " << std::endl;
            for (int j = 0; j < kruskal_tensor.size(); j++) {
                std::cout << j << "-th entry: "
                          << ((tensor_deriv1.row(j) - tensor_deriv.row(j)) / eps -
                              (tensor_hess[j] * perturb_vec).transpose())
                                 .norm()
                          << std::endl;
            }
        }
    }
}

// Test function for computing tensors per tet facet
void MiNTEnergy::TestComputeDualEdgeProjectedTensors(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd frames1 = frames;
    for (int k = 2; k <= 2; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        Eigen::VectorXd kruskal_tensor;
        Eigen::MatrixXd tensor_deriv;
        std::vector<Eigen::MatrixXd> tensor_hess;

        bool is_kruskal = true;
        int facet_id = 1;

        kruskal_tensor = ComputeBasisProjectedTensorsPerTetFace(tet_facet_dual_basis_.at(tet_id).at(facet_id), frames,
                                                                is_kruskal, tet_id, k, &tensor_deriv, &tensor_hess);

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(9);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id) = frames.row(tet_id) + eps * perturb_vec.transpose();

            Eigen::VectorXd kruskal_tensor1;
            Eigen::MatrixXd tensor_deriv1;

            // kruskal_tensor1 = ComputeKruskalTensorsPerTet(frames1, tet_id, k, &tensor_deriv1, nullptr);
            kruskal_tensor1 = ComputeBasisProjectedTensorsPerTetFace(
                tet_facet_dual_basis_.at(tet_id).at(facet_id), frames1, is_kruskal, tet_id, k, &tensor_deriv1, nullptr);

            std::cout << "\neps: " << eps << std::endl;
            std::cout << "grad-check: "
                      << ((kruskal_tensor1 - kruskal_tensor) / eps - tensor_deriv * perturb_vec).norm() << std::endl;
            std::cout << "hess-check: " << std::endl;
            for (int j = 0; j < kruskal_tensor.size(); j++) {
                std::cout << j << "-th entry: "
                          << ((tensor_deriv1.row(j) - tensor_deriv.row(j)) / eps -
                              (tensor_hess[j] * perturb_vec).transpose())
                                 .norm()
                          << std::endl;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////
///  Smoothness Energy Test Functions
///////////////////////////////////////////////////////////////////////

// Test Local Projection
void MiNTEnergy::TestPerFaceSmoothnessHessianProjection(const Eigen::MatrixXd& frames, int face_id) {
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);

    if (tet_id0 == -1 || tet_id1 == -1) {
        std::cout << "Boundary face" << std::endl;
        return;
    }

    for (int k = 2; k <= 6; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        std::vector<Eigen::VectorXd> tensors;
        std::vector<Eigen::MatrixXd> tensor_deriv;
        std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;

        tensors = ComputeKruskalTensors(frames, k, &tensor_deriv, &tensor_hess);
        Eigen::VectorXd coeff_weights = tensor_ord_weight_[k];

        Eigen::VectorXd deriv;
        Eigen::MatrixXd hess;

        double energy = ComputeSmoothnessEnergyPerFace(tensors, &tensor_deriv, &tensor_hess, coeff_weights, face_id,
                                                       &deriv, &hess, true);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(hess);
        Eigen::VectorXd D = eigensolver.eigenvalues();

        std::cout << "eigen values (expected to be >= 0):\n" << D.transpose() << std::endl;
    }
}

// Test function for computing smoothness energy per face
void MiNTEnergy::TestComputeSmoothnessEnergyPerFace(const Eigen::MatrixXd& frames, int face_id) {
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);

    if (tet_id0 == -1 || tet_id1 == -1) {
        std::cout << "Boundary face" << std::endl;
        return;
    }

    for (int k = 2; k <= 6; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        std::vector<Eigen::VectorXd> tensors;
        std::vector<Eigen::MatrixXd> tensor_deriv;
        std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;

        tensors = ComputeKruskalTensors(frames, k, &tensor_deriv, &tensor_hess);
        Eigen::VectorXd coeff_weights = tensor_ord_weight_[k];

        Eigen::VectorXd deriv;
        Eigen::MatrixXd hess;

        double energy =
            ComputeSmoothnessEnergyPerFace(tensors, &tensor_deriv, &tensor_hess, coeff_weights, face_id, &deriv, &hess);

        Eigen::MatrixXd frames1 = frames;

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(18);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id0) = frames.row(tet_id0) + eps * perturb_vec.segment<9>(0).transpose();
            frames1.row(tet_id1) = frames.row(tet_id1) + eps * perturb_vec.segment<9>(9).transpose();

            std::vector<Eigen::VectorXd> tensors1;
            std::vector<Eigen::MatrixXd> tensor1_deriv;

            tensors1 = ComputeKruskalTensors(frames1, k, &tensor1_deriv, nullptr);

            Eigen::VectorXd deriv1;
            double energy1 =
                ComputeSmoothnessEnergyPerFace(tensors1, &tensor1_deriv, nullptr, coeff_weights, face_id, &deriv1);

            std::cout << "eps: " << eps << std::endl;
            std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
            std::cout << "deriv-hess: " << ((deriv1 - deriv) / eps - hess * perturb_vec).norm() << std::endl;
        }
    }
}

// Test function for computing smoothness energy by order
void MiNTEnergy::TestComputeSmoothnessEnergyByOrder(const Eigen::MatrixXd& frames, int ord) {
    Eigen::VectorXd coeff_weights = tensor_ord_weight_[ord];

    Eigen::VectorXd deriv;
    Eigen::SparseMatrix<double> hess;

    double energy = ComputeSmoothnessEnergyByOrder(frames, ord, &deriv, &hess);

    Eigen::MatrixXd frames1 = frames;

    Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(deriv.size());

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(10, -i);

        for (int j = 0; j < mesh_->nTets(); j++) {
            frames1.row(j) = frames.row(j) + eps * perturb_vec.segment<9>(9 * j).transpose();
        }

        Eigen::VectorXd deriv1;
        double energy1 = ComputeSmoothnessEnergyByOrder(frames1, ord, &deriv1, nullptr);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
        std::cout << "deriv-hess: " << ((deriv1 - deriv) / eps - hess * perturb_vec).norm() << std::endl;
    }
}

// Test function for computing smoothness energy
void MiNTEnergy::TestComputeSmoothnessEnergy(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeSmoothnessEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting smoothness energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing smoothness energy
void MiNTEnergy::TestComputeCombedSmoothnessEnergy(const Eigen::MatrixXd& frames) {
    ComputePermutations(frames);
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeCombedSmoothnessEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting combed smoothness energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing smoothness energy
void MiNTEnergy::TestComputeCombedIntegrabilityEnergy(const Eigen::MatrixXd& frames) {
    ComputePermutations(frames);
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeCombedIntegrabilityEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting combed integrability energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing smoothness energy
void MiNTEnergy::TestComputeHodgeLaplacianSmoothnessEnergy(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeHodgeLaplacianSmoothnessEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting smoothness energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

///////////////////////////////////////////////////////////////////////
///  Integrability
///////////////////////////////////////////////////////////////////////

// Test function for computing smoothness energy per face
void MiNTEnergy::TestComputeIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id) {
    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);

    if (tet_id0 == -1 || tet_id1 == -1) {
        std::cout << "Boundary face" << std::endl;
        return;
    }

    for (int k = 2; k <= 6; k += 2) {
        std::cout << "\n========== Order " << k << " ============" << std::endl;
        std::vector<Eigen::VectorXd> tensors;
        std::vector<Eigen::MatrixXd> tensor_deriv;
        std::vector<std::vector<Eigen::MatrixXd>> tensor_hess;

        bool is_kruskal = false;

        // tensors = ComputeKruskalTensors(frames, k, &tensor_deriv, &tensor_hess);
        tensors = ComputeFaceProjectedTensors(frames, is_kruskal, k, &tensor_deriv, &tensor_hess);

        Eigen::VectorXd coeff_weights = tensor2D_ord_weight_[k];

        Eigen::VectorXd deriv;
        Eigen::MatrixXd hess;

        double energy = ComputeIntegrabilityEnergyPerFace(tensors, &tensor_deriv, &tensor_hess, coeff_weights, face_id,
                                                          &deriv, &hess);

        Eigen::MatrixXd frames1 = frames;

        Eigen::VectorXd perturb_vec = Eigen::VectorXd::Random(18);

        for (int i = 4; i < 12; i++) {
            double eps = std::pow(10, -i);
            frames1.row(tet_id0) = frames.row(tet_id0) + eps * perturb_vec.segment<9>(0).transpose();
            frames1.row(tet_id1) = frames.row(tet_id1) + eps * perturb_vec.segment<9>(9).transpose();

            std::vector<Eigen::VectorXd> tensors1;
            std::vector<Eigen::MatrixXd> tensor1_deriv;

            tensors1 = ComputeFaceProjectedTensors(frames1, is_kruskal, k, &tensor1_deriv, nullptr);

            Eigen::VectorXd deriv1;
            double energy1 =
                ComputeIntegrabilityEnergyPerFace(tensors1, &tensor1_deriv, nullptr, coeff_weights, face_id, &deriv1);

            std::cout << "eps: " << eps << std::endl;
            std::cout << "energy-deriv: " << (energy1 - energy) / eps - deriv.dot(perturb_vec) << std::endl;
            std::cout << "deriv-hess: " << ((deriv1 - deriv) / eps - hess * perturb_vec).norm() << std::endl;
        }
    }
}

void MiNTEnergy::TestIntegrabilityMintVsPrimal(const Eigen::MatrixXd& frames, bool is_kruskal) {
    double w_mint_tmp = weights_.w_mint;
    weights_.w_mint = 1e0;

    int k = 2;   // order
    // bool is_kruskal = false;
    Eigen::VectorXd coeff_weights = tensor2D_ord_weight_[k];

    std::vector<Eigen::VectorXd> tensors;
    tensors = ComputeFaceProjectedTensors(frames, is_kruskal, k, nullptr, nullptr);

    std::cout << "\nTesting Integrability Penalty" << std::endl;
    std::cout << "\nFinding facets with large discrepency between primal and symmetric integrability" << std::endl;

    // std::cout << "mint energy: " << ComputeIntegrabilityEnergyWithTriplets(frames, nullptr, nullptr, false)
    //           << std::endl;
    // std::cout << "primal energy: " << ComputePrimalIntegrabilityEnergy(frames) << std::endl;

    double largest_diff = 0;
    int largest_diff_id = 0;
    for (int curr_face = 0; curr_face < mesh_->nFaces(); curr_face++) {
        int tet_id0 = mesh_->faceTet(curr_face, 0);
        int tet_id1 = mesh_->faceTet(curr_face, 1);

        if (tet_id0 == -1 || tet_id1 == -1) {
            // std::cout << "Boundary face" << std::endl;
            continue;
        }

        double mint_energy_face_cur =
            ComputeIntegrabilityEnergyPerFace(tensors, nullptr, nullptr, coeff_weights, curr_face, nullptr, nullptr);
        double primal_energy_face_cur = ComputePrimalIntegrabilityEnergyPerFace(frames, curr_face);

        double diff = std::abs(mint_energy_face_cur - primal_energy_face_cur);

        if (diff > largest_diff) {
            largest_diff = diff;
            largest_diff_id = curr_face;
            std::cout << "face_id: " << curr_face << ", mint_energy: " << mint_energy_face_cur
                      << " primal_energy: " << primal_energy_face_cur << std::endl;
        }
    }

    int face_id = largest_diff_id;

    std::cout << "tensor computation done: " << tensors.size() << std::endl;

    int tet_id0 = mesh_->faceTet(face_id, 0);
    int tet_id1 = mesh_->faceTet(face_id, 1);

    int tetfaceidx0 = tet_id0 * 4 + mesh_->faceTetIndex(face_id, 0);
    int tetfaceidx1 = tet_id1 * 4 + mesh_->faceTetIndex(face_id, 1);

    std::cout << "tet_id0: " << tet_id0 << ", tet_id1: " << tet_id1 << std::endl;
    std::cout << "face_id: " << face_id << ", tetfaceidx0: " << tetfaceidx0 << ", tetfaceidx1: " << tetfaceidx1
              << std::endl;
    std::cout << "tensor u: " << tensors[tetfaceidx0].transpose() << ", tensor v: " << tensors[tetfaceidx1].transpose()
              << std::endl;

    Eigen::MatrixXd face_basis0 = tet_facet_basis_.at(tet_id0).at(mesh_->faceTetIndex(face_id, 0));
    Eigen::VectorXd vec0 = frames.row(tet_id0).segment(3 * 0, 3);
    Eigen::VectorXd proj_vec0 = face_basis0 * vec0;

    Eigen::VectorXd vec1 = frames.row(tet_id1).segment(3 * 0, 3);
    Eigen::VectorXd proj_vec1 = face_basis0 * vec1;

    std::cout << "full u: " << vec0.transpose() << ", full v: " << vec1.transpose() << std::endl;
    std::cout << "u: " << proj_vec0.transpose() << ", v: " << proj_vec1.transpose() << std::endl;

    double weight = dual_face_volumes_[face_id] / (dual_edge_lengths_[face_id] * dual_edge_lengths_[face_id]);
    std::cout << "primal error: "
              << (proj_vec0 - proj_vec1).squaredNorm() * (proj_vec0 - proj_vec1).squaredNorm() * weight << std::endl;

    double mint_energy_face =
        ComputeIntegrabilityEnergyPerFace(tensors, nullptr, nullptr, coeff_weights, face_id, nullptr, nullptr);
    double primal_energy_face = ComputePrimalIntegrabilityEnergyPerFace(frames, face_id);

    std::cout << "face_id: " << face_id << ", mint_energy: " << mint_energy_face
              << " primal_energy: " << primal_energy_face << std::endl;

    weights_.w_mint = w_mint_tmp;

    return;
}

// Test function for computing integrablity  energy
void MiNTEnergy::TestComputeIntegrabilityEnergy(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeIntegrabilityEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting Integrability energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing deviation energy
void MiNTEnergy::TestComputeDeviationEnergy(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeDeviationEnergy(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting deviation energy" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing unit norm penalty
void MiNTEnergy::TestComputeUnitNormalBoundaryPenalty(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeUnitNormalBoundaryPenalty(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting boundary alignment penalty" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing unit norm penalty
void MiNTEnergy::TestComputeUnitNormPenalty(const Eigen::MatrixXd& frames) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, Eigen::VectorXd* deriv,
                           Eigen::SparseMatrix<double>* hess, bool is_PSD_proj) {
        return ComputeUnitNormPenalty(input_frames, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting unit norm penalty" << std::endl;
    TestGeneralFunction(frames, energy_func);
}

// Test function for computing unit norm penalty per tet
void MiNTEnergy::TestComputeUnitNormPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int element_id, Eigen::VectorXd* deriv,
                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
        return ComputeUnitNormPenaltyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting unit norm penalty per tet" << std::endl;
    TestGeneralFunctionPerTet(frames, tet_id, energy_func);
}

// Test function for computing unit norm penalty per tet
void MiNTEnergy::TestComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd test_frames;
    test_frames = frames;
    test_frames.row(tet_id) *= 1. / (frames.row(tet_id).norm() + 1e-12);
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int element_id, Eigen::VectorXd* deriv,
                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
        return ComputeDeterminantPenaltyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting Unit Determinant penalty per tet" << std::endl;
    TestGeneralFunctionPerTet(test_frames, tet_id, energy_func);
}

// Test function for computing unit norm penalty per tet
void MiNTEnergy::TestComputeOrthogonalityEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd test_frames;
    test_frames = frames;
    test_frames.row(tet_id) *= 1. / (frames.row(tet_id).norm() + 1e-12);
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int element_id, Eigen::VectorXd* deriv,
                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
        return ComputeOrthogonalityEnergyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting Orthogonality penalty per tet" << std::endl;
    TestGeneralFunctionPerTet(test_frames, tet_id, energy_func);
}

// Test function for computing unit norm penalty per tet
void MiNTEnergy::TestComputeScaledJacobianPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd test_frames;
    test_frames = frames;
    test_frames.row(tet_id) *= 1. / (frames.row(tet_id).norm() + 1e-12);
    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int element_id, Eigen::VectorXd* deriv,
                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
        return ComputeScaledJacobianPenaltyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
        // return ComputeScaledJacobianInversePenaltyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting Scaled Jacobian penalty per tet" << std::endl;
    TestGeneralFunctionPerTet(test_frames, tet_id, energy_func);
}

// Test function for computing unit norm penalty per tet
void MiNTEnergy::TestComputeUnitBarrierPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id) {
    Eigen::MatrixXd test_frames;
    test_frames = frames;
    // test_frames.row(tet_id) *= 1. / frames.row(tet_id).norm();
    // test_frames.row(tet_id) = Eigen::VectorXd::Random(9);
    test_frames.row(tet_id) << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    auto energy_func = [&](const Eigen::MatrixXd& input_frames, int element_id, Eigen::VectorXd* deriv,
                           Eigen::MatrixXd* hess, bool is_PSD_proj) {
        // return ComputeCrossProductEnergyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);

        return ComputeUnitBarrierPenaltyPerTet(input_frames, element_id, deriv, hess, is_PSD_proj);
    };
    std::cout << "\nTesting Unit Barrier penalty per tet" << std::endl;

    std::cout << "tet_id: " << tet_id << "frame: " << test_frames.row(tet_id) << std::endl;

    TestGeneralFunctionPerTet(test_frames, tet_id, energy_func);
}
*/
}   // namespace MiNT3D