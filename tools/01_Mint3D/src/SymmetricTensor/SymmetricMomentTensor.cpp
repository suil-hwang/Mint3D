#include "SymmetricMomentTensor.h"
#include <array>
#include <iostream>

namespace SymmetricTensor {
// Get the coefficients of the k-th order tensor of v
Eigen::VectorXd GetMomentCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv,
                                      std::vector<Eigen::MatrixXd>* hess) {
    switch (ord) {
        case 2:
            return Get2ndOrderMomentCoefficients(v, deriv, hess);
        case 4:
            return Get4thOrderMomentCoefficients(v, deriv, hess);
        case 6:
            return Get6thOrderMomentCoefficients(v, deriv, hess);
        default:
            throw std::runtime_error("The order of tensor product is not supported.");
    }
}

// Get the coefficients of the k-th order tensor of normalized v
Eigen::VectorXd GetNormalizedMomentCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv,
                                                std::vector<Eigen::MatrixXd>* hess) {
    Eigen::MatrixXd v_hat_deriv;
    std::vector<Eigen::MatrixXd> v_hat_hess;
    Eigen::VectorXd v_hat =
        GetNormalizedVector(v, (deriv || hess) ? &v_hat_deriv : nullptr, hess ? &v_hat_hess : nullptr);

    Eigen::MatrixXd deriv_hat;
    std::vector<Eigen::MatrixXd> hess_hat;

    Eigen::VectorXd coeff_hat =
        GetMomentCoefficients(v_hat, ord, (deriv || hess) ? &deriv_hat : nullptr, hess ? &hess_hat : nullptr);

    if (deriv) {
        *deriv = deriv_hat * v_hat_deriv;
    }

    if (hess) {
        int dim = v.size();
        hess->resize(coeff_hat.size());
        for (int i = 0; i < coeff_hat.size(); i++) {
            (*hess)[i] = v_hat_deriv.transpose() * hess_hat[i] * v_hat_deriv;

            for (int l = 0; l < dim; l++) {
                (*hess)[i] += deriv_hat(i, l) * v_hat_hess[l];
            }
        }
    }

    return coeff_hat;
}

//////////// 2nd order tensor ////////////
// Get the coefficients of the 2nd order tensor of v.
// In 2d, the coefficients are {v_0^2, v_0v_1, v_1^2}
// In 3d, the coefficients are {v_0^2, v_0v_1, v_0v_2, v_1^2, v_1v_2, v_2^2}
Eigen::VectorXd Get2ndOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv,
                                              std::vector<Eigen::MatrixXd>* hess) {
    int dim = v.size();
    int nvars = GetCoeffSize(dim, 2);

    Eigen::VectorXd coeff(nvars);
    if (deriv) {
        deriv->setZero(nvars, dim);
    }

    if (hess) {
        hess->resize(nvars, Eigen::MatrixXd::Zero(dim, dim));
    }

    int idx = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            coeff(idx) = v(i) * v(j);

            if (deriv) {
                (*deriv)(idx, i) += v(j);
                (*deriv)(idx, j) += v(i);
            }

            if (hess) {
                (*hess)[idx](i, j) += 1;
                (*hess)[idx](j, i) += 1;
            }

            idx++;
        }
    }
    return coeff;
}

//////////// 4th order tensor ////////////
// Get the coefficients of the 4th order tensor of v.
// In 2d, the coefficients are {v_0^4, v_0^3v_1, v_0^2v_1^2, v_0v_1^3, v_1^4}
// In 3d, the coefficients are {v_0^4, v_0^3v_1, v_0^3v_2, v_0^2v_1^2, v_0^2v_1v_2, v_0^2v_2^2, v_0v_1^3, v_0v_1^2v_2,
// v_0v_1v_2^2, v_0v_2^3, v_1^4, v_1^3v_2, v_1^2v_2^2, v_1v_2^3, v_2^4}
Eigen::VectorXd Get4thOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv,
                                              std::vector<Eigen::MatrixXd>* hess) {
    int dim = v.size();
    int nvars = GetCoeffSize(dim, 4);

    Eigen::VectorXd coeff(nvars);

    if (deriv) {
        deriv->setZero(nvars, dim);
    }

    if (hess) {
        hess->resize(nvars, Eigen::MatrixXd::Zero(dim, dim));
    }

    int idx = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            for (int m = j; m < dim; m++) {
                for (int n = m; n < dim; n++) {
                    coeff(idx) = v(i) * v(j) * v(m) * v(n);

                    if (deriv) {
                        (*deriv)(idx, i) += v(j) * v(m) * v(n);
                        (*deriv)(idx, j) += v(i) * v(m) * v(n);
                        (*deriv)(idx, m) += v(i) * v(j) * v(n);
                        (*deriv)(idx, n) += v(i) * v(j) * v(m);
                    }

                    if (hess) {
                        (*hess)[idx](i, j) += v(m) * v(n);
                        (*hess)[idx](j, i) += v(m) * v(n);

                        (*hess)[idx](i, m) += v(j) * v(n);
                        (*hess)[idx](m, i) += v(j) * v(n);

                        (*hess)[idx](i, n) += v(j) * v(m);
                        (*hess)[idx](n, i) += v(j) * v(m);

                        (*hess)[idx](j, m) += v(i) * v(n);
                        (*hess)[idx](m, j) += v(i) * v(n);

                        (*hess)[idx](j, n) += v(i) * v(m);
                        (*hess)[idx](n, j) += v(i) * v(m);

                        (*hess)[idx](m, n) += v(i) * v(j);
                        (*hess)[idx](n, m) += v(i) * v(j);
                    }

                    idx++;
                }
            }
        }
    }

    return coeff;
}

//////////// 6th order tensor ////////////
// Get the coefficients of the 6th order tensor of v.
// In 2d, the coefficients are {v_0^6, v_0^5v_1, v_0^4v_1^2, v_0^3v_1^3, v_0^2v_1^4, v_0v_1^5, v_1^6}
// In 3d, the coefficients are {v_0^6, v_0^5v_1, v_0^5v_2, v_0^4v_1^2, v_0^4v_1v_2, v_0^4v_2^2, v_0^3v_1^3,
// v_0^3v_1^2v_2, v_0^3v_1v_2^2, v_0^3v_2^3, v_0^2v_1^4, v_0^2v_1^3v_2, v_0^2v_1^2v_2^2, v_0^2v_1v_2^3, v_0^2v_2^4,
// v_0v_1^5, v_0v_1^4v_2, v_0v_1^3v_2^2, v_0v_1^2v_2^3, v_0v_1v_2^4, v_0v_2^5, v_1^6, v_1^5v_2, v_1^4v_2^2, v_1^3v_2^3,
// v_1^2v_2^4, v_1v_2^5, v_2^6}
Eigen::VectorXd Get6thOrderMomentCoefficients(const Eigen::VectorXd& v, Eigen::MatrixXd* deriv,
                                              std::vector<Eigen::MatrixXd>* hess) {
    int dim = v.size();
    int nvars = GetCoeffSize(dim, 6);

    Eigen::VectorXd coeff(nvars);

    if (deriv) {
        deriv->setZero(nvars, dim);
    }

    if (hess) {
        hess->resize(nvars, Eigen::MatrixXd::Zero(dim, dim));
    }

    int idx = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            for (int k = j; k < dim; k++) {
                for (int l = k; l < dim; l++) {
                    for (int m = l; m < dim; m++) {
                        for (int n = m; n < dim; n++) {
                            coeff(idx) = v(i) * v(j) * v(k) * v(l) * v(m) * v(n);
                            if (deriv) {
                                (*deriv)(idx, i) += v(j) * v(k) * v(l) * v(m) * v(n);
                                (*deriv)(idx, j) += v(i) * v(k) * v(l) * v(m) * v(n);
                                (*deriv)(idx, k) += v(i) * v(j) * v(l) * v(m) * v(n);
                                (*deriv)(idx, l) += v(i) * v(j) * v(k) * v(m) * v(n);
                                (*deriv)(idx, m) += v(i) * v(j) * v(k) * v(l) * v(n);
                                (*deriv)(idx, n) += v(i) * v(j) * v(k) * v(l) * v(m);
                            }
                            if (hess) {
                                (*hess)[idx](i, j) += v(k) * v(l) * v(m) * v(n);
                                (*hess)[idx](j, i) += v(k) * v(l) * v(m) * v(n);

                                (*hess)[idx](i, k) += v(j) * v(l) * v(m) * v(n);
                                (*hess)[idx](k, i) += v(j) * v(l) * v(m) * v(n);

                                (*hess)[idx](i, l) += v(j) * v(k) * v(m) * v(n);
                                (*hess)[idx](l, i) += v(j) * v(k) * v(m) * v(n);

                                (*hess)[idx](i, m) += v(j) * v(k) * v(l) * v(n);
                                (*hess)[idx](m, i) += v(j) * v(k) * v(l) * v(n);

                                (*hess)[idx](i, n) += v(j) * v(k) * v(l) * v(m);
                                (*hess)[idx](n, i) += v(j) * v(k) * v(l) * v(m);

                                (*hess)[idx](j, k) += v(i) * v(l) * v(m) * v(n);
                                (*hess)[idx](k, j) += v(i) * v(l) * v(m) * v(n);

                                (*hess)[idx](j, l) += v(i) * v(k) * v(m) * v(n);
                                (*hess)[idx](l, j) += v(i) * v(k) * v(m) * v(n);

                                (*hess)[idx](j, m) += v(i) * v(k) * v(l) * v(n);
                                (*hess)[idx](m, j) += v(i) * v(k) * v(l) * v(n);

                                (*hess)[idx](j, n) += v(i) * v(k) * v(l) * v(m);
                                (*hess)[idx](n, j) += v(i) * v(k) * v(l) * v(m);

                                (*hess)[idx](k, l) += v(i) * v(j) * v(m) * v(n);
                                (*hess)[idx](l, k) += v(i) * v(j) * v(m) * v(n);

                                (*hess)[idx](k, m) += v(i) * v(j) * v(l) * v(n);
                                (*hess)[idx](m, k) += v(i) * v(j) * v(l) * v(n);

                                (*hess)[idx](k, n) += v(i) * v(j) * v(l) * v(m);
                                (*hess)[idx](n, k) += v(i) * v(j) * v(l) * v(m);

                                (*hess)[idx](l, m) += v(i) * v(j) * v(k) * v(n);
                                (*hess)[idx](m, l) += v(i) * v(j) * v(k) * v(n);

                                (*hess)[idx](l, n) += v(i) * v(j) * v(k) * v(m);
                                (*hess)[idx](n, l) += v(i) * v(j) * v(k) * v(m);

                                (*hess)[idx](m, n) += v(i) * v(j) * v(k) * v(l);
                                (*hess)[idx](n, m) += v(i) * v(j) * v(k) * v(l);
                            }
                            idx++;
                        }
                    }
                }
            }
        }
    }
    return coeff;
}

// Test Coeffient derivatives and hessians
void TestMomentCoefficientDerivs(const Eigen::VectorXd& v, int ord) {
    Eigen::VectorXd rand_v = Eigen::VectorXd::Random(v.size());

    Eigen::MatrixXd deriv;
    std::vector<Eigen::MatrixXd> hess;
    Eigen::VectorXd coeff = GetMomentCoefficients(v, ord, &deriv, &hess);

    std::cout << "\n=== order: " << ord << ", vec: " << v.transpose() << std::endl;
    std::cout << "coeff " << coeff.transpose() << std::endl;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(0.1, i);
        Eigen::MatrixXd deriv1;
        Eigen::VectorXd v1 = v + eps * rand_v;
        Eigen::VectorXd coeff1 = GetMomentCoefficients(v1, ord, &deriv1);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "grad-check: " << ((coeff1 - coeff) / eps - deriv * rand_v).norm() << std::endl;
        std::cout << "hess-check: " << std::endl;
        for (int j = 0; j < coeff.size(); j++) {
            std::cout << j
                      << "-th entry: " << ((deriv1.row(j) - deriv.row(j)) / eps - (hess[j] * rand_v).transpose()).norm()
                      << std::endl;
        }
    }
}

// Test Coeffient derivatives and hessians of normalized moment tensor
void TestNormalizedMomentCoefficientDerivs(const Eigen::VectorXd& v, int ord) {
    Eigen::VectorXd rand_v = Eigen::VectorXd::Random(v.size());

    Eigen::MatrixXd deriv;
    std::vector<Eigen::MatrixXd> hess;
    Eigen::VectorXd coeff = GetNormalizedMomentCoefficients(v, ord, &deriv, &hess);

    std::cout << "\n=== order: " << ord << ", vec: " << v.transpose() << std::endl;
    std::cout << "coeff " << coeff.transpose() << std::endl;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(0.1, i);
        Eigen::MatrixXd deriv1;
        Eigen::VectorXd v1 = v + eps * rand_v;
        Eigen::VectorXd coeff1 = GetNormalizedMomentCoefficients(v1, ord, &deriv1);

        std::cout << "\neps: " << eps << std::endl;
        std::cout << "grad-check: " << ((coeff1 - coeff) / eps - deriv * rand_v).norm() << std::endl;
        std::cout << "hess-check: " << std::endl;
        for (int j = 0; j < coeff.size(); j++) {
            std::cout << j
                      << "-th entry: " << ((deriv1.row(j) - deriv.row(j)) / eps - (hess[j] * rand_v).transpose()).norm()
                      << std::endl;
        }
    }
}
}   // namespace SymmetricTensor
