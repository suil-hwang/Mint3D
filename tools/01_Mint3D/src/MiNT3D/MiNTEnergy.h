#pragma once

#include <Eigen/Dense>
/*
 * The class to compute the MiNT energy of a frame field on a tetrahedral mesh
 * Author: Zhen Chen
 * There are several things to optimize:
 * 1. It might be better to write a unified function to compute the energy, and then call it for different energy types
 * 2. We can write a unified function to assemmble the hessian from per-tet hessians
 */

#include <Eigen/Sparse>

#include "CubeCover/FrameField.h"
#include "CubeCover/TetMeshConnectivity.h"
#include "Surface.h"

#include "MiNTSolveParams.h"

#include <iomanip>   // Required for std::setprecision and std::scientific
#include <iostream>

namespace MiNT3D {

class MiNTEnergy {
   public:
    MiNTEnergy() {
        V_ = nullptr;
        mesh_ = nullptr;
        mesh_not_extended_ = nullptr;
        dual_edge_lengths_.resize(0);
        dual_face_volumes_ = {};
        tensor_ord_weight_.clear();
        tensor2D_ord_weight_.clear();
        tet_volumes_ = {};
        show_energy_ = false;
    }
    ~MiNTEnergy() {}

    // The constructor, where V is the vertex matrix and mesh is the tetrahedral mesh connectivity
    MiNTEnergy(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
               const CubeCover::TetMeshConnectivity& mesh_orig)
        : V_(&V), mesh_(&mesh), mesh_not_extended_(&mesh_orig) {
        permutation_as_smooth_as_possible_.resize(mesh_->nFaces());
        permutation_as_integrable_as_possible_.resize(mesh_->nFaces());

        show_energy_ = false;
    }

    // Precompute the data structures for the energy computation, including the dual volumes of the faces, and dual
    // vertex positions (currently we use the volume centroid as the dual vertex position)
    void Precompute();

    // Compute an orthogonal basis per tet facet
    void ComputePerTetFacetBasis();

    // Compute an orthogonal basis per tet facet
    void ComputeBoundaryNormals();

    void ComputePermutations(const Eigen::MatrixXd& frames);

    // Compute Guided Tensors
    void ComputeGuidedTensors(const Eigen::MatrixXd& guided_frames);

    /// smoothness
    // Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron
    double ComputeSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                   Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    double ComputeCombedSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                         Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    double ComputeCombedIntegrabilityEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                            Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                               std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                               bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeCombedSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                     std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                     bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeCombedIntegrabilityEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                        std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                        bool is_PSD_proj = false, bool is_kruskal = false);

    /// hodge laplacian smoothness
    // Compute the smoothness energy of the frame field by minimizing |div f|^2 + |curl f|^2
    // . We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron
    double ComputeHodgeLaplacianSmoothnessEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                 Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the smoothness energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeHodgeLaplacianSmoothnessEnergyWithTriplets(
        const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
        std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr, bool is_PSD_proj = false);

    /// integrability
    // Compute the integrability energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron
    double ComputeIntegrabilityEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                      Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false,
                                      bool is_kruskal = false);

    // Compute the integrability energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeIntegrabilityEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                  std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                  bool is_PSD_proj = false, bool is_kruskal = false,
                                                  double weight = -1.);

    // Compute the primal integrability energy per face.
    // the MiNT solver should succeed in driving this to a small value if the integrability constraints are on and
    // everything is working correctly.
    double ComputePrimalIntegrabilityEnergy(const Eigen::MatrixXd& frames);

    /// divergence
    // Compute the divergence energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron
    double ComputeDivergenceEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                   Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false,
                                   bool is_kruskal = false);

    // Compute the divergence energy of the frame field. We assume that frames.rows() == mesh_->nTets() and
    // frames.cols() == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeDivergenceEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                               std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                               bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the primal divergence energy per face.
    double ComputePrimalDivergenceEnergy(const Eigen::MatrixXd& frames);

    /// per tet energies

    // Fitting term
    // Compute the deviation energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols()
    // == 9, that is 3 frames per tetrahedron
    double ComputeDeviationEnergy(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                  Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the deviation energy of the frame field. We assume that frames.rows() == mesh_->nTets() and frames.cols()
    // == 9, that is 3 frames per tetrahedron, but we return the hessian in triplet form
    double ComputeDeviationEnergyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                              std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                              bool is_PSD_proj = false);

    // Unit Norm
    // Compute the unit norm penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
    // frames per tetrahedron
    double ComputeUnitNormPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                  Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeUnitNormPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                              std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                              bool is_PSD_proj = false);

    // Soft boundary alignment
    // Compute the unit normal boundary penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9,
    // that is 3 frames per tetrahedron
    double ComputeUnitNormalBoundaryPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                            Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Soft boundary condition
    double ComputeUnitNormalBoundaryPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                        std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                        bool is_PSD_proj = false);

    // Scaled Jacobian
    // Compute the unit norm penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
    // frames per tetrahedron
    double ComputeScaledJacobianPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                        Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeScaledJacobianPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                    std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                    bool is_PSD_proj = false);

    // Viscosity
    // Compute the unit norm penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
    // frames per tetrahedron
    double ComputeViscosityPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                   Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeViscosityPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                               std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                               bool is_PSD_proj = false);

    // Feature Alignment
    // add a penalty to the ghost tets to align with sharp features.
    double ComputeFeatureAlignmentPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                          Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeFeatureAlignmentPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                      std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                      bool is_PSD_proj = false);

    // Unit Determinant Penalty
    // Compute the unit determinant penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that
    // is 3 frames per tetrahedron
    double ComputeDeterminantPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                     Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeDeterminantPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                 std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                 bool is_PSD_proj = false);

    // Unit Norm Barrier
    // Compute the unit norm penalty. We assume that frames.rows() == mesh_->nTets() and frames.cols() == 9, that is 3
    // frames per tetrahedron
    double ComputeUnitBarrierPenalty(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                     Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty with the triplet
    double ComputeUnitBarrierPenaltyWithTriplets(const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv = nullptr,
                                                 std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                 bool is_PSD_proj = false);

    ////////////////////// Unit tests //////////////////////
    // Test function for general function
    void TestGeneralFunction(
        const Eigen::MatrixXd& frames,
        std::function<double(const Eigen::MatrixXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)>
            energy_func);

    // Test function for general function defined per tet
    void TestGeneralFunctionPerTet(
        const Eigen::MatrixXd& frames, int tet_id,
        std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func);

    // Test function for general function defined per face
    void TestGeneralFunctionPerFace(
        const Eigen::MatrixXd& frames, int face_id,
        std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func);

    /// helpers
    // Test Local Projection
    void TestPerFaceSmoothnessHessianProjection(const Eigen::MatrixXd& frames, int face_id);

    // Test function for computing Kruskal tensors per tet
    void TestComputeMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // Test function for computing Kruskal tensors per tet
    void TestComputeKruskalTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // Test function for computing projected tensors per face
    void TestComputeFaceProjectedTensors(const Eigen::MatrixXd& frames, int tet_id);

    // Test function for computing projected tensors per dual edge
    void TestComputeDualEdgeProjectedTensors(const Eigen::MatrixXd& frames, int tet_id);

    /// smoothness
    // Test function for computing smoothness energy
    void TestComputeSmoothnessEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing smoothness energy
    void TestComputeCombedSmoothnessEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing smoothness energy
    void TestComputeCombedIntegrabilityEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing smoothness energy
    void TestComputeHodgeLaplacianSmoothnessEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing smoothness energy by order
    void TestComputeSmoothnessEnergyByOrder(const Eigen::MatrixXd& frames, int ord);

    // Test function for computing smoothness energy per face
    void TestComputeSmoothnessEnergyPerFace(const Eigen::MatrixXd& frames, int face_id);

    /// integrability
    // Test function for computing integrability energy
    void TestComputeIntegrabilityEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing integrability energy by order
    void TestComputeIntegrabilityEnergyByOrder(const Eigen::MatrixXd& frames, int ord);

    // Test function for computing integrability energy per face
    void TestComputeIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id);

    void TestIntegrabilityMintVsPrimal(const Eigen::MatrixXd& frames, bool is_kruskal);

    /// Divergence
    // Test function for computing divergence energy
    void TestComputeDivergenceEnergy(const Eigen::MatrixXd& frames);

    // Test function for computing divergence energy by order
    void TestComputeDivergenceEnergyByOrder(const Eigen::MatrixXd& frames, int ord);

    // Test function for computing divergence energy per face
    void TestComputeDivergenceEnergyPerFace(const Eigen::MatrixXd& frames, int face_id);

    void TestDivergenceMintVsPrimal(const Eigen::MatrixXd& frames, bool is_kruskal);

    /// per tet energies

    // Test function for computing deviation energy
    void TestComputeDeviationEnergy(const Eigen::MatrixXd& frames);

    // Test Unit Boundary Penalty
    void TestComputeUnitNormalBoundaryPenalty(const Eigen::MatrixXd& frames);

    // Test Unit norm energy
    void TestComputeUnitNormPenalty(const Eigen::MatrixXd& frames);
    // Test unit norm energy per tet
    void TestComputeUnitNormPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // Test Scaled Jacobian energy
    void TestScaledJacobianNormPenalty(const Eigen::MatrixXd& frames);
    // Test scaled jacobian energy per tet
    void TestComputeScaledJacobianPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // Test scaled jacobian energy per tet
    void TestComputeViscosityPenalty(const Eigen::MatrixXd& frames);

    // Test feature alignment energy per tet
    void TestComputeFeatureAlignmentPenalty(const Eigen::MatrixXd& frames);

    // Test orthogonality energy per tet
    void TestComputeOrthogonalityEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // // Test Scaled Jacobian energy
    // void TestDeterminantNormPenalty(const Eigen::MatrixXd& frames);
    // Test scaled jacobian energy per tet
    void TestComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id);

    // Test Unit Barrier energy
    void TestUnitBarrierPenalty(const Eigen::MatrixXd& frames);
    // Test Unit Barrier energy per tet
    void TestComputeUnitBarrierPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id);

   public:
    // Assembling to the total energy from per tet computation
    double AssembleEnergyFromTets(
        const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv, std::vector<Eigen::Triplet<double>>* hess_triplets,
        bool is_PSD_proj,
        std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func);

    // Assembling to the total energy from per face computation
    double AssembleEnergyFromFaces(
        const Eigen::MatrixXd& frames, Eigen::VectorXd* deriv, std::vector<Eigen::Triplet<double>>* hess_triplets,
        bool is_PSD_proj,
        std::function<double(const Eigen::MatrixXd&, int, Eigen::VectorXd*, Eigen::MatrixXd*, bool)> energy_func);

    // Compute the Kruskal tensors
    std::vector<Eigen::VectorXd> ComputeMomentTensors(const Eigen::MatrixXd& frames, int order,
                                                      std::vector<Eigen::MatrixXd>* deriv = nullptr,
                                                      std::vector<std::vector<Eigen::MatrixXd>>* hess = nullptr);

    // Compute the Moment tensor of a frame field on a tetrahedron
    Eigen::VectorXd ComputeMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                               Eigen::MatrixXd* deriv = nullptr,
                                               std::vector<Eigen::MatrixXd>* hess = nullptr);

    // Compute the Kruskal tensors
    std::vector<Eigen::VectorXd> ComputeKruskalTensors(const Eigen::MatrixXd& frames, int order,
                                                       std::vector<Eigen::MatrixXd>* deriv = nullptr,
                                                       std::vector<std::vector<Eigen::MatrixXd>>* hess = nullptr);

    // Compute the Kruskal tensor of a frame field on a tetrahedron
    Eigen::VectorXd ComputeKruskalTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                Eigen::MatrixXd* deriv = nullptr,
                                                std::vector<Eigen::MatrixXd>* hess = nullptr);

    // Compute the frames projected to each tet face and then lifted as krushkal tensors
    std::vector<Eigen::VectorXd> ComputeFaceProjectedTensors(const Eigen::MatrixXd& frames, bool is_kruskal, int order,
                                                             std::vector<Eigen::MatrixXd>* deriv = nullptr,
                                                             std::vector<std::vector<Eigen::MatrixXd>>* hess = nullptr);

    // Compute the frames projected to each tet face and then lifted as krushkal tensors
    std::vector<Eigen::VectorXd> ComputeDualEdgeProjectedTensors(
        const Eigen::MatrixXd& frames, bool is_kruskal, int order, std::vector<Eigen::MatrixXd>* deriv = nullptr,
        std::vector<std::vector<Eigen::MatrixXd>>* hess = nullptr);

    // Compute the Kruskal tensor of a frame field projected to a tet face
    Eigen::VectorXd ComputeBasisProjectedTensorsPerTetFace(const Eigen::MatrixXd& tet_facet_basis,
                                                           const Eigen::MatrixXd& frames, bool is_kruskal, int tet_id,
                                                           int order, Eigen::MatrixXd* deriv = nullptr,
                                                           std::vector<Eigen::MatrixXd>* hess = nullptr);

    // Compute the Hodge Laplacian energy of the frame field by order.
    /*
    * Input:
    * frames: the current frames
    * ord: the order of the Hodge Laplacian
    * deriv: the derivative of the Hodge Laplacian energy in terms of tet frames
    * hess: the hessian of the Hodge Laplacian energy in terms of tet frames
    * is_PSD_proj: whether to project the hessian to the PSD cone
    *
    * Output:
    * the Hodge Laplacian energy
    *
    * Requirement:
    * The precomputed data structures should be initialized by Precompute()
    *
    * Return:
    * the Hodge Laplacian energy
    *
    * Note:
    * The Hodge Laplacian energy is computed as the finite difference of Kruskal tensors from the two adjacent tets,
    * weighted by the dual volume of the face.
    *



    */
    // Input:

    double ComputeHodgeLaplacianEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                          Eigen::VectorXd* deriv,
                                                          std::vector<Eigen::Triplet<double>>* hess_triplets,
                                                          bool is_PSD_proj);

    // Compute the combed smoothness energy of the frame field per face
    // For each face, the smoothness energy is computed as the finite difference of combed vectors from the two
    // adjacent tets, weighted by the dual volume of the face.
    /*
     * Input:
     * frames: the current frame variables
     * permutations: per edge permutation of the frame field
     * face_id: the face id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the smoothness energy in terms of tet frames
     * hess: the hessian of the smoothness energy in terms of tet frames
     *
     * Return:
     * the smoothness energy per face
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeCombedSmoothnessEnergyPerFace(const Eigen::MatrixXd& frames, int face_id,
                                                Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                bool is_PSD_proj = false);

    // Compute the combed integrability energy of the frame field per face
    // For each face, the integrability energy is computed as the finite difference of combed vectors from the two
    // adjacent tets, weighted by the dual volume of the face.
    /*
     * Input:
     * frames: the current frame variables
     * permutations: per edge permutation of the frame field
     * face_id: the face id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the integrability energy in terms of tet frames
     * hess: the hessian of the integrability energy in terms of tet frames
     *
     * Return:
     * the integrability energy per face
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeCombedIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id,
                                                   Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                   bool is_PSD_proj = false);

    // Defined above
    // // Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice
    // that
    // // we always append to the triplets instead of clearing it.
    // double ComputeCombedSmoothnessEnergyWithTriplets(const Eigen::MatrixXd& frames, int ord,
    //                                                   Eigen::VectorXd* deriv = nullptr,
    //                                                   std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
    //                                                   bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the smoothness energy of the frame field per face by order
    // For each face, the smoothness energy is computed as the finite difference of Kruskal tensors from the two
    // adjacent tets, weighted by the dual volume of the face.
    /*
     * Input:
     * tensors: the Kruskal tensors of the frame field on the tets
     * tensor_deriv: the derivatives of the Kruskal tensors in terms of tet frames
     * tensor_hess: the hessians of the Kruskal tensors in terms of tet frames
     * coeff_weights: the weight of the Kruskal tensor, which is used to massure the norm of the tensor difference
     * face_id: the face id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the smoothness energy in terms of tet frames
     * hess: the hessian of the smoothness energy in terms of tet frames
     *
     * Return:
     * the smoothness energy per face
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeSmoothnessEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                          std::vector<Eigen::MatrixXd>* tensor_deriv,
                                          std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                          const Eigen::VectorXd& coeff_weights, int face_id,
                                          Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                          bool is_PSD_proj = false);

    double ComputeScaleFreeSmoothnessEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                                   std::vector<Eigen::MatrixXd>* tensor_deriv,
                                                   std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                                   const Eigen::VectorXd& coeff_weights, int face_id, double p,
                                                   Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                   bool is_PSD_proj = false);

    // Compute the smoothness energy of the frame field by order.
    double ComputeSmoothnessEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv = nullptr,
                                          Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false,
                                          bool is_kruskal = false);

    // Compute the smoothness energy of the frame field by order, but we return the hessian in triplet form. Notice that
    // we always append to the triplets instead of clearing it.
    double ComputeSmoothnessEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                      Eigen::VectorXd* deriv = nullptr,
                                                      std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                      bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the integrability energy of the frame field per face by order
    // For each face, the integrability energy is computed as the finite difference of moment tensors from the two
    // adjacent tets, projected onto the dual facet and weighted by the dual volume of the face.
    /*
     * Input:
     * tensors: the moment tensors of the frame field on the tets
     * tensor_deriv: the derivatives of the moment tensors in terms of tet frames
     * tensor_hess: the hessians of the moment tensors in terms of tet frames
     * coeff_weights: the weight of the moment tensor, which is used to measure the norm of the tensor difference
     * face_id: the face id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the smoothness energy in terms of tet frames
     * hess: the hessian of the smoothness energy in terms of tet frames
     *
     * Return:
     * the symmetric integrability residual per face
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeIntegrabilityEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                             std::vector<Eigen::MatrixXd>* tensor_deriv,
                                             std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                             const Eigen::VectorXd& coeff_weights, int face_id,
                                             Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                             bool is_PSD_proj = false);

    // Compute the integrability energy of the frame field by order.
    double ComputeIntegrabilityEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv = nullptr,
                                             Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false,
                                             bool is_kruskal = false);

    // Compute the integrability energy of the frame field by order, but we return the hessian in triplet form. Notice
    // that we always append to the triplets instead of clearing it.
    double ComputeIntegrabilityEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                         Eigen::VectorXd* deriv = nullptr,
                                                         std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                         bool is_PSD_proj = false, bool is_kruskal = false,
                                                         double weight = -1.);

    // Compute the primal integrability energy per face.
    // the MiNT solver should succeed in driving this to a small value
    // if the integrability constraints are on and everything is working correctly.
    double ComputePrimalIntegrabilityEnergyPerFace(const Eigen::MatrixXd& frames, int face_id);

    // Compute the divergence energy of the frame field per face by order
    // For each face, the divergence energy is computed as the finite difference of moment tensors from the two
    // adjacent tets, projected onto the dual facet and weighted by the dual volume of the face.
    /*
     * Input:
     * tensors: the moment tensors of the frame field on the tets
     * tensor_deriv: the derivatives of the moment tensors in terms of tet frames
     * tensor_hess: the hessians of the moment tensors in terms of tet frames
     * coeff_weights: the weight of the moment tensor, which is used to measure the norm of the tensor difference
     * face_id: the face id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the smoothness energy in terms of tet frames
     * hess: the hessian of the smoothness energy in terms of tet frames
     *
     * Return:
     * the symmetric divergence residual per face
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeDivergenceEnergyPerFace(const std::vector<Eigen::VectorXd>& tensors,
                                          std::vector<Eigen::MatrixXd>* tensor_deriv,
                                          std::vector<std::vector<Eigen::MatrixXd>>* tensor_hess,
                                          const Eigen::VectorXd& coeff_weights, int face_id,
                                          Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                          bool is_PSD_proj = false);

    // Compute the Divergence energy of the frame field by order.
    double ComputeDivergenceEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv = nullptr,
                                          Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false,
                                          bool is_kruskal = false);

    // Compute the Divergence energy of the frame field by order, but we return the hessian in triplet form. Notice
    // that we always append to the triplets instead of clearing it.
    double ComputeDivergenceEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                      Eigen::VectorXd* deriv = nullptr,
                                                      std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                      bool is_PSD_proj = false, bool is_kruskal = false);

    // Compute the primal Divergence energy per face.
    // the MiNT solver should succeed in driving this to a small value
    // if the Divergence constraints are on and everything is working correctly.
    double ComputePrimalDivergenceEnergyPerFace(const Eigen::MatrixXd& frames, int face_id);

    // Compute the normalized moment tensors
    std::vector<Eigen::VectorXd> ComputeNormalizedMomentTensors(
        const Eigen::MatrixXd& frames, int order, std::vector<Eigen::MatrixXd>* deriv = nullptr,
        std::vector<std::vector<Eigen::MatrixXd>>* hess = nullptr);

    // Compute the normalized moment tensor of a frame field on a tetrahedron
    Eigen::VectorXd ComputeNormalizedMomentTensorsPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                                         Eigen::MatrixXd* deriv = nullptr,
                                                         std::vector<Eigen::MatrixXd>* hess = nullptr);

    // Compute the deviation energy per tet by order
    /*
     * Input:
     * frames: the current frames
     * face_id: the tet id
     * order: the order of the normalized moment tensor
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the deviation energy in terms of tet frames
     * hess: the hessian of the deviation energy in terms of tet frames
     *
     * Return:
     * the deviation energy per tet
     *
     * Requirement:
     * The precomputed data structures should be initialized by Precompute()
     */
    double ComputeDeviationEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, int order,
                                        Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                        bool is_PSD_proj = false);

    // Compute the deviation energy of the frame field by order.
    double ComputeDeviationEnergyByOrder(const Eigen::MatrixXd& frames, int ord, Eigen::VectorXd* deriv = nullptr,
                                         Eigen::SparseMatrix<double>* hess = nullptr, bool is_PSD_proj = false);

    // Compute the deviation energy of the frame field by order, but we return the hessian in triplet form. Notice that
    // we always append to the triplets instead of clearing it.
    double ComputeDeviationEnergyByOrderWithTriplets(const Eigen::MatrixXd& frames, int ord,
                                                     Eigen::VectorXd* deriv = nullptr,
                                                     std::vector<Eigen::Triplet<double>>* hess_triplets = nullptr,
                                                     bool is_PSD_proj = false);

    double ComputeOrthogonalityEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                            Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

    // Compute the unit norm penalty per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the unit norm penalty in terms of tet frames
     * hess: the hessian of the unit norm penalty in terms of tet frames
     *
     * Return:
     * the unit norm penalty per tet
     */
    double ComputeUnitNormPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                        Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

    double ComputeUnitNormalBoundaryPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                                  Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                  bool is_PSD_proj = false);

    //
    double ComputeDeterminantPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                           Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

    // Compute the scaled jacobian penalty per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the scaled jacobian penalty in terms of tet frames
     * hess: the hessian of the scaled jacobian penalty in terms of tet frames
     *
     * Return:
     * the scaled jacobian penalty per tet
     */
    double ComputeScaledJacobianPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                              Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                              bool is_PSD_proj = false);

    double ComputeScaledJacobianInversePenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                                     Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                     bool is_PSD_proj = false);

    // Compute the viscosity penalty per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * context (not passed in), the state to enforce L2 loss too.
     * Hacky, sorry, lazy.
     *
     * Output:
     * deriv: the derivative of the viscosity penalty in terms of tet frames
     * hess: the hessian of the viscosity penalty in terms of tet frames
     *
     * Return:
     * the viscosity penalty per tet
     */
    double ComputeViscosityPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                         Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

    // Compute the feature alignment penalty per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * context (not passed in), the state to enforce L2 loss too.
     * Hacky, sorry, lazy.
     *
     * Output:
     * deriv: the derivative of the feature alignment penalty in terms of tet frames
     * hess: the hessian of the feature alignment penalty in terms of tet frames
     *
     * Return:
     * the feature alignment penalty per tet
     */
    double ComputeFeatureAlignmentPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id,
                                                std::vector<Eigen::Vector3d>& tet_features,
                                                Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr,
                                                bool is_PSD_proj = false);

    // Compute the unit barrier penalty per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the scaled jacobian penalty in terms of tet frames
     * hess: the hessian of the scaled jacobian penalty in terms of tet frames
     *
     * Return:
     * the penalty of (1 - 1 / (min_i |f_i| ) )^2
     */
    double ComputeUnitBarrierPenaltyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                           Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

    //           ComputeCrossProductEnergyPerTet

    // Compute the cross product energy per tet
    /*
     * Input:
     * frames: the current frames
     * tet_id: the tet id
     * is_PSD_proj: whether to project the hessian to the PSD cone
     *
     * Output:
     * deriv: the derivative of the scaled jacobian penalty in terms of tet frames
     * hess: the hessian of the scaled jacobian penalty in terms of tet frames
     *
     * Return:
     * the penalty of (1 - 1 / (min_i |f_i| ) )^2
     */

    double ComputeCrossProductEnergyPerTet(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv = nullptr,
                                           Eigen::MatrixXd* hess = nullptr, bool is_PSD_proj = false);

   public:
    MiNT3D::SolveParams weights_;
    std::unordered_map<int, Eigen::VectorXd> tensor_ord_weight_;     // the weight of the Kruskal tensor of each order
    std::unordered_map<int, Eigen::VectorXd> tensor2D_ord_weight_;   // the weight of the Kruskal tensor of each order

   public:
    const Eigen::MatrixXd* V_;                                  // the vertex position matrix pointer
    const CubeCover::TetMeshConnectivity* mesh_;                // the tetrahedral mesh connectivity pointer
    const CubeCover::TetMeshConnectivity* mesh_not_extended_;   // the tetrahedral mesh connectivity pointer

    // Choice of combing operation might change how these are computed and updated.
    // This formulation was introduced to avoid the potential issue
    // where multiple symmetric terms do not pull back to a consistent
    // low rank frame field.
    // Instead we fix the combing of all operators non-symmetric operators consistently.
    // for now, this is set when you call combed smoothness.

    std::vector<Eigen::MatrixXd> permutation_as_smooth_as_possible_;
    std::vector<Eigen::MatrixXd> permutation_as_integrable_as_possible_;

    bool show_energy_;
    Eigen::VectorXd per_tet_energy_density_;

    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_basis_;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_dual_basis_;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_full_basis_;

    std::vector<std::vector<Eigen::Vector3d>> boundary_features_;

    Eigen::MatrixXd boundary_normals_;
    Eigen::MatrixXd boundary_b1_;
    Eigen::MatrixXd boundary_b2_;

    Eigen::MatrixXd frames_prev_outer_step_;
    Eigen::MatrixXd frames_prev_inner_step_;

    std::vector<double> dual_edge_lengths_;   // the dual lengths of the edges
    std::vector<double> dual_face_volumes_;   // the dual volumes of the faces
    std::vector<double> tet_volumes_;         // the volumes of the tets
    double ave_edge_len_;   // the average edge length, this is used to make all the energy has the unit of length

    std::unordered_map<int, std::vector<Eigen::VectorXd>>
        guided_tensor_coeffs_;   // the weight of the guided Kruskal tensor of each order

    std::unordered_map<int, std::vector<Eigen::VectorXd>>
        self_align_tensor_coeffs_;   // the weight of the guided Kruskal tensor of each order
};
}   // namespace MiNT3D