// #include "Projection.h"

// namespace MiNT3D {


// // constraint mask 
// Projection::Projection(const std::vector<bool>& keep_dofs) {
//     int fulldofs = keep_dofs.size();
//     dofmap.resize(fulldofs);
//     int idx = 0;
//     constraint_mask_ = Eigen::VectorXd::Zero(fulldofs);
//     for (int i = 0; i < fulldofs; i++) {
//         if (keep_dofs[i]) {
//             dofmap[i] = idx;
//             invdofmap.push_back(i);
//             idx++;
//         } else {
//             dofmap[i] = -1;
//             constraint_mask_[i] = 1;
//         }
//     }

    
// }

// void Projection::setFixedVars(const Eigen::MatrixXd& ext_frames) {
//     // multiply coefficient wise with the constraint mask
    
//     fixed_vars_ = MatrixToVector(ext_frames);

//     fixed_vars_ = fixed_vars_.cwiseProduct(constraint_mask_);
// }

// int Projection::ProjDOFs() const {
//     return invdofmap.size();
// }

// void Projection::ProjectVector(const Eigen::VectorXd& full_vec, Eigen::VectorXd& proj_vec) const {
//     int projdofs = invdofmap.size();
//     proj_vec.resize(projdofs);

//     for (int i = 0; i < projdofs; i++) {
//         proj_vec[i] = full_vec[invdofmap[i]];
//     }
// }

// void Projection::UnprojectVector(const Eigen::VectorXd& proj_vec, Eigen::VectorXd& full_vec) const {
//     int fulldofs = dofmap.size();
//     full_vec.resize(fulldofs);
//     for (int i = 0; i < fulldofs; i++) {
//         full_vec[i] = (dofmap[i] == -1 ? 0.0 : proj_vec[dofmap[i]]);
//     }
// }

// void Projection::ProjectMatrix(std::vector<Eigen::Triplet<double> >& mat) const {
//     int dim = mat.size();
//     for (int i = 0; i < dim; i++) {
//         int r = mat[i].row();
//         int c = mat[i].col();
//         int pr = dofmap[r];
//         int pc = dofmap[c];
//         if (pr != -1 && pc != -1) {
//             mat[i] = {pr, pc, mat[i].value()};
//         } else {
//             mat[i] = {0, 0, 0.0};
//         }
//     }
// }


// // Convert extended frames to variable (moved from MiNTModel)
// Eigen::VectorXd Projection::ConvertFramesToVariable(const Eigen::MatrixXd& ext_frames) const {
//     Eigen::VectorXd full_x = ext_frames.reshaped();  // Reshape to a vector
//     Eigen::VectorXd proj_x;
//     this->ProjectVector(full_x, proj_x);
//     return proj_x;
// }

// // Convert variable to frames (moved from MiNTModel)
// Eigen::MatrixXd Projection::ConvertVariableToFrames(const Eigen::VectorXd& x) const {
//     Eigen::VectorXd full_x;
//     this->UnprojectVector(x, full_x);
//     full_x += fixed_vars_;

//     // Ensure correct size for reshaping
//     assert(full_x.size() % 9 == 0);  

//     int num_frames = full_x.size() / 9;
//     return Eigen::Map<Eigen::MatrixXd>(full_x.data(), num_frames, 9); // Reshape to a matrix
// }

// // Set the selection matrix (moved from MiNTModel)
// Eigen::SparseMatrix<double> Projection::ComputeSelectionMatrix(const std::vector<bool>& fixed_dofs) {
//     std::vector<Eigen::Triplet<double>> T;
//     int full_dofs = fixed_dofs.size();
//     int rows = 0;
//     for (int i = 0; i < full_dofs; i++) {
//         if (!fixed_dofs[i]) { // Free variable
//             T.push_back(Eigen::Triplet<double>(rows, i, 1));
//             rows++;
//         }
//     }
//     Eigen::SparseMatrix<double> select_mat(rows, full_dofs);
//     select_mat.setFromTriplets(T.begin(), T.end());
//     return select_mat;
// }


// // Flat matrix to vector (moved from MiNTModel)
// Eigen::VectorXd Projection::MatrixToVector(const Eigen::MatrixXd& mat) {
//     Eigen::VectorXd vec(mat.rows() * mat.cols());
//     for (int i = 0; i < mat.rows(); i++) {
//         vec.segment(i * mat.cols(), mat.cols()) = mat.row(i);
//     }
//     return vec;
// }

// // Vector to matrix (moved from MiNTModel)
// Eigen::MatrixXd Projection::VectorToMatrix(const Eigen::VectorXd& vec, int cols) {
//     Eigen::MatrixXd mat(vec.size() / cols, cols);
//     for (int i = 0; i < vec.size() / cols; i++) {
//         mat.row(i) = vec.segment(i * cols, cols);
//     }
//     return mat;
// }


// } // namespace MiNT3D
