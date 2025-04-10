#pragma once

#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "CubeCover/TetMeshConnectivity.h"

#include "../Optimization/NewtonDescent.h"
#include "MiNTEnergy.h"
#include "MiNTSolveParams.h"

namespace MiNT3D {
class Projection {
   public:
    Projection() : is_porj_needed_(true) {}
    Projection(const std::vector<bool>& keep_dofs);

    // the number of dofs after projection
    int ProjDOFs() const { return invdofmap.size(); }

    // project the full vector to the reduced vector
    void ProjectVector(const Eigen::VectorXd& full_vec, Eigen::VectorXd& proj_vec) const;

    // unproject the reduced vector to the full vector with fixed dofs set to 0
    void UnprojectVector(const Eigen::VectorXd& proj_vec, Eigen::VectorXd& full_vec) const;

    // project the full matrix to the reduced matrix
    void ProjectMatrix(std::vector<Eigen::Triplet<double>>& mat) const;

    void setFixedVars(const Eigen::MatrixXd& ext_frames);
    Eigen::VectorXd constraint_mask_;
    Eigen::VectorXd fixed_vars_;

   private:
    std::vector<int> dofmap;
    std::vector<int> invdofmap;
    bool is_porj_needed_;
};

class EnergyLog {
   public:
    std::vector<double> smoothness;
    std::vector<double> unit_penalty;
    std::vector<double> symmetric_integrability;
    std::vector<double> primal_integrability;
    std::vector<double> scaled_jacobian;
    std::vector<double> combed_integrability;
    std::vector<double> asap_combed_smoothness;
    std::vector<double> aiap_combed_smoothness;

    std::vector<double> combed_smoothness;
    std::vector<double> unit_barrier;
    std::vector<double> boundary_penalty;

    std::vector<double> curl_weight;
    std::vector<double> total_energy_unscale;

    std::vector<double> fit_penalty;

    std::vector<double> global_scale_outiter;

    bool logToFile(std::string logPath);
};

class MiNTModel {
   public:
    MiNTModel() {
        init_V_ = nullptr;
        init_mesh_ = nullptr;
        extended_V_ = Eigen::MatrixXd(0, 3);
        extended_mesh_ = CubeCover::TetMeshConnectivity();
        lambda_integrability_ = 0.0;
        lambda_regularity_ = 0.0;
        lambda_scalar_jacobi_ = 0.0;

        proj_ = nullptr;

        solveData = std::unique_ptr<OptSolver::SolveData>(new OptSolver::SolveData());   // new SolveData();
        energyLog = std::unique_ptr<EnergyLog>(new EnergyLog());                         // new SolveData();
    }

    ~MiNTModel() = default;

    // Set the Mesh
    void SetMesh(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh) {
        init_V_ = &V;
        init_mesh_ = &mesh;
    }

    void SetSharpFeatureEdges(const std::vector<std::vector<Eigen::Vector3d>>& sharp_feature_edges) {
        sharp_feature_edges_ = sharp_feature_edges;
    }

    // Convert variable to frames, defined on the extended mesh
    Eigen::MatrixXd ConvertVariableToFrames(const Eigen::VectorXd& x);

    // Convert the extended frames to variable
    Eigen::VectorXd ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames);

    // Evaluate the energy
    double ComputeEnergy(const Eigen::MatrixXd& input_frames, Eigen::MatrixXd* virtual_ext_bnd_frames = nullptr);

    // Solve the smooth frames along with boundary normals
    bool SolveSmoothFrames(SolveParams solve_params, const Eigen::MatrixXd& input_frames,
                           Eigen::MatrixXd& output_frames, Eigen::MatrixXd* virtual_ext_bnd_frames = nullptr,
                           Eigen::MatrixXd* opt_virtual_ext_bnd_frames = nullptr,
                           std::function<void()>* viz_callback = nullptr);

    // hack to show the density of a particular energy.
    bool ShowDensity(SolveParams solve_params, const Eigen::MatrixXd& input_frames, Eigen::MatrixXd& output_frames,
                     Eigen::MatrixXd* virtual_ext_bnd_frames, Eigen::MatrixXd* opt_virtual_ext_bnd_frames);

    bool LogStepState(const Eigen::VectorXd& x, const SolveParams& solve_params, std::string logPath);

    bool LogSolveData(const OptSolver::SolveData* solveData, const EnergyLog* energyLog, std::string logPath);

    std::unique_ptr<OptSolver::SolveData> solveData;
    std::unique_ptr<EnergyLog> energyLog;

   public:
    // Set the extended mesh. Basically, the extended mesh is the mesh with additional layers for the boundary by
    // flipping to properly handle the boundary conditions
    void AddAdditionalLayerForBoundary();

    // Set the selection matrix
    Eigen::SparseMatrix<double> ComputeSelectionMatrix(const std::vector<bool>& fixed_dofs);

    // Flat matrix to vector
    static Eigen::VectorXd MatrixToVector(const Eigen::MatrixXd& mat);

    // Vector to matrix
    static Eigen::MatrixXd VectorToMatrix(const Eigen::VectorXd& vec, int cols);

   public:
    Eigen::MatrixXd cur_smooth_density;
    Eigen::MatrixXd cur_mint_density;
    Eigen::MatrixXd cur_neohookean_density;
    Eigen::MatrixXd cur_asap_smoothness_density;
    Eigen::MatrixXd cur_aiap_smoothness_density;

    Eigen::MatrixXd cur_unit_norm_density;
    Eigen::MatrixXd cur_bound_align_density;
    Eigen::MatrixXd cur_scaled_jacobian_density;

    Eigen::VectorXd flagged_frames;
    Eigen::VectorXd flagged_edges;

    Eigen::MatrixXd cur_primal_integrability_density;

   public:
    const Eigen::MatrixXd* init_V_;                     // the initial vertex matrix
    const CubeCover::TetMeshConnectivity* init_mesh_;   // the initial tetrahedral mesh connectivity
    Eigen::MatrixXd extended_V_;                        // the extended vertex matrix
    Eigen::MatrixXi
        extended_T_;   // the extended tetrahedral cell matrix <-- only used for init extended mesh and visualization
    CubeCover::TetMeshConnectivity extended_mesh_;   // the tetrahedral mesh connectivity

    Eigen::VectorXd fixed_normal_frames_;     // the fixed normal frames
    Eigen::VectorXd fixed_interior_frames_;   // the fixed normal frames

    std::unique_ptr<Projection> proj_;   // the projection object

    std::vector<std::vector<Eigen::Vector3d>> sharp_feature_edges_;

    int currentFileID_ = 0;
    std::string logInnerIter;
    std::string logOuterIter;

    // std::string exp_name = "";
    std::string logPath;

    double lambda_integrability_;   // the penalty parameter for the integrability constraint
    double lambda_regularity_;      // the penalty parameter for the regularity constraint
    double lambda_scalar_jacobi_;   // the penalty parameter for the scalar Jacobi
};
}   // namespace MiNT3D