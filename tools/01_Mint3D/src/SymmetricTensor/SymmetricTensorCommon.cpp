#include "SymmetricTensorCommon.h"
#include <iostream>

namespace SymmetricTensor {
double GetNorm(const Eigen::VectorXd& v, Eigen::VectorXd* deriv, Eigen::MatrixXd* hess) {
    double norm = v.norm();
    if (deriv) {
        *deriv = v / norm;
    }
    if (hess) {
        *hess = Eigen::MatrixXd::Zero(v.size(), v.size());
        for (int i = 0; i < v.size(); i++) {
            for (int j = 0; j < v.size(); j++) {
                if (i == j) {
                    (*hess)(i, j) = (norm * norm - v(i) * v(i)) / (norm * norm * norm);
                } else {
                    (*hess)(i, j) = -v(i) * v(j) / (norm * norm * norm);
                }
            }
        }
    }
    return norm;
}

// Get the normalized vector
Eigen::VectorXd GetNormalizedVector(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv,
                                    std::vector<Eigen::MatrixXd>* hess) {
    int dim = v.size();
    Eigen::VectorXd v_hat;
    double norm = GetNorm(v, &v_hat, deriv);

    if (hess) {
        hess->resize(dim, Eigen::MatrixXd::Zero(dim, dim));
        double norm3 = norm * norm * norm;
        double norm5 = norm3 * norm * norm;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    (*hess)[i](j, k) = 3 * v(i) * v(j) * v(k) / norm5;
                    if (j == k) {
                        (*hess)[i](j, k) -= v(i) / norm3;
                    }
                    if (k == i) {
                        (*hess)[i](j, k) -= v(j) / norm3;
                    }
                    if (i == j) {
                        (*hess)[i](j, k) -= v(k) / norm3;
                    }
                }
            }
        }
    }

    return v_hat;
}

// Get the weight of coefficients
Eigen::VectorXd GetCoeffWeights(int dim, int ord) {
    switch (ord) {
        case 2:
            return Get2ndOrderCoeffWeights(dim);
        case 4:
            return Get4thOrderCoeffWeights(dim);
        case 6:
            return Get6thOrderCoeffWeights(dim);
        default:
            throw std::runtime_error("The order of tensor product is not supported.");
    }
}

// Get the size of the coefficient vector. It is dim + ord - 1 choose ord.
int GetCoeffSize(int dim, int ord) {
    int numerate = 1;
    int denominate = 1;

    for (int i = 0; i < ord; i++) {
        numerate *= dim + i;
        denominate *= i + 1;
    }

    return numerate / denominate;
}

// Get the weight of coefficients of the 2nd order tensor.
// In 2d, the weights are {1, 2, 1}
// In 3d, the weights are {1, 2, 2, 1, 2, 1}
Eigen::VectorXd Get2ndOrderCoeffWeights(int dim) {
    Eigen::VectorXd weights;
    switch (dim) {
        case 1: {
            weights.resize(1);
            weights << 1;
            break;
        }
        case 2: {
            weights.resize(3);
            weights << 1, 2, 1;
            break;
        }
        case 3: {
            weights.resize(6);
            weights << 1, 2, 2, 1, 2, 1;
            // weights *= 1e-1;
            // weights *= 1e-2;
            break;
        }
        default:
            throw std::runtime_error("The dimension is not supported.");
    }
    return weights;
}

// Get the weight of coefficients of the 4th order tensor.
// In 2d, the weights are {1, 4, 6, 4, 1}
// In 3d, the weights are {1, 4, 4, 6, 12, 6, 4, 12, 12, 4, 1, 4, 6, 4, 1}
Eigen::VectorXd Get4thOrderCoeffWeights(int dim) {
    Eigen::VectorXd weights;
    switch (dim) {
        case 2: {
            weights.resize(5);
            weights << 1, 4, 6, 4, 1;
            break;
        }
        case 3: {
            weights.resize(15);
            weights << 1, 4, 4, 6, 12, 6, 4, 12, 12, 4, 1, 4, 6, 4, 1;
            // weights *= 1e-1; 
            // weights *= .5;
            break;
        }
        default:
            throw std::runtime_error("The dimension is not supported.");
    }
    return weights;
}

// Get the weight of coefficients of the 6th order tensor.
// In 2d, the weights are {1, 6, 15, 20, 15, 6, 1}
// In 3d, the weights are {1, 6, 6, 15, 30, 15, 20, 60, 60, 20, 15, 60, 90, 60, 15, 6, 30, 60, 60, 30, 6, 1, 6, 15, 20,
// 15, 6, 1}
Eigen::VectorXd Get6thOrderCoeffWeights(int dim) {
    Eigen::VectorXd weights;
    switch (dim) {
        case 2: {
            weights.resize(7);
            weights << 1, 6, 15, 20, 15, 6, 1;
            break;
        }
        case 3: {
            weights.resize(28);
            weights << 1, 6, 6, 15, 30, 15, 20, 60, 60, 20, 15, 60, 90, 60, 15, 6, 30, 60, 60, 30, 6, 1, 6, 15, 20, 15,
                6, 1;
            // weights *= 1e1;

            // weights *= 10000.;
            break;
        }
        default:
            throw std::runtime_error("The dimension is not supported.");
    }
    return weights;
}

// Get n!
int Factorial(int n) {
    int fact = 1;
    for (int i = 2; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

// Test norm derivatives and hessians
void TestNormDerivs(const Eigen::VectorXd& v) {
    Eigen::VectorXd rand_v = Eigen::VectorXd::Random(v.size());

    Eigen::VectorXd deriv;
    Eigen::MatrixXd hess;
    double norm = GetNorm(v, &deriv, &hess);

    std::cout << "\n=== vec: " << v.transpose() << std::endl;
    std::cout << "|v| " << norm << std::endl;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(0.1, i);
        Eigen::VectorXd deriv1;
        Eigen::VectorXd v1 = v + eps * rand_v;
        double norm1 = GetNorm(v1, &deriv1);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "grad-check: " << ((norm1 - norm) / eps - deriv.dot(rand_v)) << std::endl;
        std::cout << "hess-check: " << std::endl;
        std::cout << ((deriv1 - deriv) / eps - hess * rand_v).norm() << std::endl;
    }
}

// Test normalized vector derivatives and hessians
void TestNormalizedVectorDerivs(const Eigen::VectorXd& v) {
    Eigen::VectorXd rand_v = Eigen::VectorXd::Random(v.size());

    Eigen::MatrixXd deriv;
    std::vector<Eigen::MatrixXd> hess;
    Eigen::VectorXd v_hat = GetNormalizedVector(v, &deriv, &hess);

    std::cout << "\n=== vec: " << v.transpose() << std::endl;
    std::cout << "v_hat " << v_hat.transpose() << std::endl;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(0.1, i);
        Eigen::MatrixXd deriv1;
        Eigen::VectorXd v1 = v + eps * rand_v;
        Eigen::VectorXd v_hat1 = GetNormalizedVector(v1, &deriv1);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "grad-check: " << ((v_hat1 - v_hat) / eps - deriv * rand_v).norm() << std::endl;
        std::cout << "hess-check: " << std::endl;
        for (int j = 0; j < v.size(); j++) {
            std::cout << j
                      << "-th entry: " << ((deriv1.row(j) - deriv.row(j)) / eps - (hess[j] * rand_v).transpose()).norm()
                      << std::endl;
        }
    }
}

// Test whether the weights are correct
void TestWeights(int dim, int ord) {
    switch (ord) {
        case 2: {
            Test2ndWeights(dim);
            break;
        }
        case 4: {
            Test4thWeights(dim);
            break;
        }
        case 6: {
            Test6thWeights(dim);
            break;
        }
    }
}

// Test 2nd order weights
void Test2ndWeights(int dim) {
    std::cout << "\n=== Test hardcoded weight for dim: " << dim << std::endl;
    Eigen::VectorXd weights = Get2ndOrderCoeffWeights(dim);
    int nvars = GetCoeffSize(dim, 2);
    Eigen::VectorXd weights_ref(nvars);

    int idx = 0;
    int fact = Factorial(2);
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            std::vector<int> count(dim, 0);
            count[i]++;
            count[j]++;

            weights_ref(idx) = fact;
            for (int p = 0; p < dim; p++) {
                weights_ref(idx) /= Factorial(count[p]);
            }
            idx++;
        }
    }
    std::cout << "2nd order weights: " << weights.transpose() << std::endl;
    std::cout << "2nd order weights ref: " << weights_ref.transpose() << std::endl;
    std::cout << "2nd order weights diff: " << (weights - weights_ref).norm() << std::endl;
}

// Test 4th order weights
void Test4thWeights(int dim) {
    std::cout << "\n=== Test hardcoded weight for dim: " << dim << std::endl;
    Eigen::VectorXd weights = Get4thOrderCoeffWeights(dim);
    int nvars = GetCoeffSize(dim, 4);
    Eigen::VectorXd weights_ref(nvars);

    int idx = 0;
    int fact = Factorial(4);
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            for (int m = j; m < dim; m++) {
                for (int n = m; n < dim; n++) {
                    std::vector<int> count(dim, 0);
                    count[i]++;
                    count[j]++;
                    count[m]++;
                    count[n]++;

                    weights_ref(idx) = fact;
                    for (int p = 0; p < dim; p++) {
                        weights_ref(idx) /= Factorial(count[p]);
                    }
                    idx++;
                }
            }
        }
    }

    std::cout << "4th order weights: " << weights.transpose() << std::endl;
    std::cout << "4th order weights ref: " << weights_ref.transpose() << std::endl;
    std::cout << "4th order weights diff: " << (weights - weights_ref).norm() << std::endl;
}

// Test 6th order weights
void Test6thWeights(int dim) {
    std::cout << "\n=== Test hardcoded weight for dim: " << dim << std::endl;
    Eigen::VectorXd weights = Get6thOrderCoeffWeights(dim);
    int nvars = GetCoeffSize(dim, 6);
    Eigen::VectorXd weights_ref(nvars);

    int idx = 0;
    int fact = Factorial(6);
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            for (int k = j; k < dim; k++) {
                for (int l = k; l < dim; l++) {
                    for (int m = l; m < dim; m++) {
                        for (int n = m; n < dim; n++) {
                            std::vector<int> count(dim, 0);
                            count[i]++;
                            count[j]++;
                            count[k]++;
                            count[l]++;
                            count[m]++;
                            count[n]++;

                            weights_ref(idx) = fact;
                            for (int p = 0; p < dim; p++) {
                                weights_ref(idx) /= Factorial(count[p]);
                            }
                            idx++;
                        }
                    }
                }
            }
        }
    }
    std::cout << "6th order weights: " << weights.transpose() << std::endl;
    std::cout << "6th order weights ref: " << weights_ref.transpose() << std::endl;
    std::cout << "6th order weights diff: " << (weights - weights_ref).norm() << std::endl;
}
}   // namespace SymmetricTensor