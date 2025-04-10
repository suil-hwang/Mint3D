#pragma once
#include "SymmetricTensorCommon.h"

/*
 * The Kruskal tensor is a symmetric tensor that can be represented as the sum of the tensor product of vectors.
 * We provide the Kruskal tensor of 2nd, 4th, and 6th order, where the dimension of the vector is 2 or 3, where
 * the Kruskal tensor := ||v|| \hat{v}^{\otimes ord}. (v is a vector, \hat{v} is the unit vector of v, and \otimes is
 * the tensor product) Author: Zhen Chen
 */

namespace SymmetricTensor {
// Get the coefficients of the k-th order symmetric moment tensor of v
Eigen::VectorXd GetKruskalCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv = nullptr,
                                       std::vector<Eigen::MatrixXd>* hess = nullptr);

// Test Coeffient derivatives and hessians
void TestKruskalCoefficientDerivs(const Eigen::VectorXd& v, int ord);

};   // namespace SymmetricTensor