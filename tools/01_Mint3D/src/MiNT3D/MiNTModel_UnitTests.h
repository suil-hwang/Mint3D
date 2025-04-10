#pragma once

#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "CubeCover/TetMeshConnectivity.h"

#include "../Optimization/NewtonDescent.h"
#include "MiNTEnergy.h"
#include "MiNTModel.h"

#include "MiNTSolveParams.h"

namespace MiNT3D {

class MiNTModel_UnitTests {
   public:
    MiNTModel_UnitTests() {}

    ~MiNTModel_UnitTests() = default;

    // Set the Mesh
    void SetMesh(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh) {
        init_V_ = &V;
        init_mesh_ = &mesh;
    }

    // Convert variable to frames, defined on the extended mesh
    Eigen::MatrixXd ConvertVariableToFrames(const Eigen::VectorXd& x);

    // Convert the extended frames to variable
    Eigen::VectorXd ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames);

    // Solve the smooth frames along with boundary normals
    bool SolveTwoTetTest(MiNT3D::SolveParams solve_params, const Eigen::MatrixXd& input_frames,
                         Eigen::MatrixXd& output_frames);

    std::unique_ptr<OptSolver::SolveData> solveData;
    std::unique_ptr<EnergyLog> energyLog;

   private:
    // Flat matrix to vector
    Eigen::VectorXd MatrixToVector(const Eigen::MatrixXd& mat);

    // Vector to matrix
    Eigen::MatrixXd VectorToMatrix(const Eigen::VectorXd& vec, int cols);

   public:
    const Eigen::MatrixXd* init_V_;                     // the initial vertex matrix
    const CubeCover::TetMeshConnectivity* init_mesh_;   // the initial tetrahedral mesh connectivity

    Eigen::VectorXd fixed_normal_frames_;   // the fixed normal frames
    std::unique_ptr<Projection> proj_;      // the projection object
};
}   // namespace MiNT3D
