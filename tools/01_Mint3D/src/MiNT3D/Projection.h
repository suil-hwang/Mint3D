// #ifndef PROJECTION_H
// #define PROJECTION_H

// #include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <vector>

// namespace MiNT3D {

// class Projection {
// public:
//     Projection(const std::vector<bool>& keep_dofs);

//     int ProjDOFs() const;
//     void ProjectVector(const Eigen::VectorXd& full_vec, Eigen::VectorXd& proj_vec) const;
//     void UnprojectVector(const Eigen::VectorXd& proj_vec, Eigen::VectorXd& full_vec) const;
//     void ProjectMatrix(std::vector<Eigen::Triplet<double> >& mat) const;
//     Eigen::VectorXd ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames) const;
//     Eigen::MatrixXd ConvertVariableToFrames(const Eigen::VectorXd& x) const;

//     // Methods for matrix/vector conversions
//     static Eigen::SparseMatrix<double> ComputeSelectionMatrix(const std::vector<bool>& fixed_dofs);
//     static Eigen::VectorXd MatrixToVector(const Eigen::MatrixXd& mat);
//     static Eigen::MatrixXd VectorToMatrix(const Eigen::VectorXd& vec, int cols);

//     void setFixedVars(const Eigen::MatrixXd& ext_frames);
//     Eigen::VectorXd constraint_mask_;
//     Eigen::VectorXd fixed_vars_;

    
// private:
//     std::vector<int> dofmap;
//     std::vector<int> invdofmap;
// };

// } // namespace MiNT3D

// #endif // PROJECTION_H
