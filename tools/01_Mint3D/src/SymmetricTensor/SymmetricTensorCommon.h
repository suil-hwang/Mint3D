#pragma once
#include <Eigen/Dense>
#include <vector>

namespace SymmetricTensor {
// Get the norm of a vector
double GetNorm(const Eigen::VectorXd& v, Eigen::VectorXd* deriv = nullptr, Eigen::MatrixXd* hess = nullptr);

// Test norm derivatives and hessians
void TestNormDerivs(const Eigen::VectorXd& v);

// Get the normalized vector
Eigen::VectorXd GetNormalizedVector(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv = nullptr,
                                    std::vector<Eigen::MatrixXd>* hess = nullptr);

// Test normalized vector derivatives and hessians
void TestNormalizedVectorDerivs(const Eigen::VectorXd& v);

// Get the weight of coefficients
Eigen::VectorXd GetCoeffWeights(int dim, int ord);

// Get the size of the coefficient vector
int GetCoeffSize(int dim, int ord);

// Get n!
int Factorial(int n);

// Test whether the weights are correct
void TestWeights(int dim, int ord);
// Test 2nd order weights
void Test2ndWeights(int dim);
// Test 4th order weights
void Test4thWeights(int dim);
// Test 6th order weights
void Test6thWeights(int dim);

//////////// 2nd order tensor ////////////
// Get the weight of coefficients of the 2nd order tensor
Eigen::VectorXd Get2ndOrderCoeffWeights(int dim);

//////////// 4th order tensor ////////////
// Get the weight of coefficients of the 4th order tensor
Eigen::VectorXd Get4thOrderCoeffWeights(int dim);

//////////// 6th order tensor ////////////
// Get the weight of coefficients of the 6th order tensor
Eigen::VectorXd Get6thOrderCoeffWeights(int dim);
}   // namespace SymmetricTensor