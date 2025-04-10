#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

namespace MiNT3D {
// Self-adjoint SPD projection for a small size matrix A
Eigen::MatrixXd SelfAdjSPDProjection(const Eigen::MatrixXd& A);

// Get determinant of a 3x3 matrix A, togiether with its derivative and hessian in terms of A00, A01, ..., A22
double Det2x3(const Eigen::Matrix3d& A, Eigen::VectorXd* deriv = nullptr, Eigen::Matrix<double, 9, 9>* hess = nullptr);

// Get determinant of a 3x3 matrix A, togiether with its derivative and hessian in terms of A00, A01, ..., A22
double Det3x3(const Eigen::Matrix3d& A, Eigen::VectorXd* deriv = nullptr, Eigen::Matrix<double, 9, 9>* hess = nullptr);

// Helper function to create a cross product matrix
Eigen::Matrix3d cross_product_matrix(const Eigen::Vector3d& v);

Eigen::MatrixXd ComputeOptimalPermutationMatrix(const Eigen::VectorXd& f, const Eigen::VectorXd& g,
                                                const Eigen::MatrixXd& basis);
void Test_ComputeOptimalPermutationMatrix();

Eigen::MatrixXd ExpandFramePermutation(const Eigen::MatrixXd& perm, const Eigen::MatrixXd& sgn);

Eigen::MatrixXd tensorProduct(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B);

// Test the derivative and hessian of Det3x3
void TestDet3x3();
}   // namespace MiNT3D