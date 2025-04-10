#pragma once
#include "SymmetricTensorCommon.h"

/*
 * This is the class to perform tensor product of a single vector v. We provide
 * v\otimes v, v\otimes v\otimes v\otimes v, and v^{\otimes 6}, namely, 2nd, 4th, and 6th order tensors.
 * We treat these tensors as the coefficients of the tensor product of the basis, where we always assume that the basis
 * are othonormal. For instance, in 3D, let v = \sum c_i b_i, then v\otimes v = \sum c_i c_j b_i\otimes b_j, where
 * b_i\otimes b_j is the tensor product of the basis. then v\otimes v is equivalent to the vector {c_1^2, c_1c_2,
 * c_1c_3, c_2c_1, c_2^2, c_2c_3, c_3c_1, c_3c_2, c_3^2} in the basis b_i\otimes b_j. Notice that c_1 c_2 = c_2 c_1, so
 * the order of the coefficients does not matter, although the order of the basis does matter. {c_1^2, c_1c_2, c_1c_3,
 * c_2c_1, c_2^2, c_2c_3, c_3c_1, c_3c_2, c_3^2} is equivalent to a 3x3 symmetric matrix, and we only store the unique
 * coefficients in a vector as {c_1^2, c_1c_2, c_1c_3, c_2^2, c_2c_3, c_3^2}, together with the number of counts of each
 * unique coefficient as {1, 2, 2, 1, 2, 1}. (Consider (v1 + v2 + v3)^2)
 *
 * Currently, I only provide the 2nd, 4th, and 6th order tensors, where dim(v) = 2 or 3, but it is easy to extend to
 * higher order tensors. Author: Zhen Chen
 */

namespace SymmetricTensor {
// Get the coefficients of the k-th order symmetric moment tensor of v
Eigen::VectorXd GetMomentCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv = nullptr,
                                      std::vector<Eigen::MatrixXd>* hess = nullptr);

// Get the coefficients of the k-th order symmetric moment tensor of normalized v
Eigen::VectorXd GetNormalizedMomentCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv = nullptr,
                                                std::vector<Eigen::MatrixXd>* hess = nullptr);

// Test Coeffient derivatives and hessians
void TestMomentCoefficientDerivs(const Eigen::VectorXd& v, int ord);

// Test Coeffient derivatives and hessians of normalized moment tensor
void TestNormalizedMomentCoefficientDerivs(const Eigen::VectorXd& v, int ord);

//////////// 2nd order tensor ////////////
// Get the coefficients of the 2nd order tensor of v
Eigen::VectorXd Get2ndOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv = nullptr,
                                              std::vector<Eigen::MatrixXd>* hess = nullptr);

//////////// 4th order tensor ////////////
// Get the coefficients of the 4th order tensor of v
Eigen::VectorXd Get4thOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv = nullptr,
                                              std::vector<Eigen::MatrixXd>* hess = nullptr);

//////////// 6th order tensor ////////////
// Get the coefficients of the 6th order tensor of v
Eigen::VectorXd Get6thOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv = nullptr,
                                              std::vector<Eigen::MatrixXd>* hess = nullptr);

};   // namespace SymmetricTensor