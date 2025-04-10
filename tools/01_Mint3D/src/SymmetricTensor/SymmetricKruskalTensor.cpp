#include "SymmetricKruskalTensor.h"

#include <array>
#include <iostream>
#include <vector>

#include <tbb/parallel_for.h>

#include "SymmetricMomentTensor.h"

namespace SymmetricTensor {
// Get the coefficients of the k-th order tensor of v
Eigen::VectorXd GetKruskalCoefficients(const Eigen::VectorXd& v, int ord, Eigen::MatrixXd* deriv,
                                       std::vector<Eigen::MatrixXd>* hess) {
    Eigen::MatrixXd v_hat_deriv;
    std::vector<Eigen::MatrixXd> v_hat_hess;
    Eigen::VectorXd v_hat =
        GetNormalizedVector(v, (deriv || hess) ? &v_hat_deriv : nullptr, hess ? &v_hat_hess : nullptr);

    Eigen::VectorXd norm_deriv;
    Eigen::MatrixXd norm_hess;
    double norm = GetNorm(v, (deriv || hess) ? &norm_deriv : nullptr, hess ? &norm_hess : nullptr);

    Eigen::MatrixXd deriv_hat;
    std::vector<Eigen::MatrixXd> hess_hat;

    Eigen::VectorXd coeff_hat =
        GetMomentCoefficients(v_hat, ord, (deriv || hess) ? &deriv_hat : nullptr, hess ? &hess_hat : nullptr);

    Eigen::VectorXd coeff = norm * coeff_hat;

    if (deriv) {
        *deriv = coeff_hat * norm_deriv.transpose() + norm * deriv_hat * v_hat_deriv;
    }

    if (hess) {
        int dim = v.size();
        hess->resize(coeff.size());
        //            tbb::parallel_for(
        //                tbb::blocked_range<int>(0, coeff.size()),
        //                [&](const tbb::blocked_range<int>& r) {
        //                    for (int i = r.begin(); i < r.end(); i++) {
        //                        Eigen::MatrixXd tmp_mat = norm_deriv * deriv_hat.row(i) * v_hat_deriv;
        //                        (*hess)[i] = coeff_hat[i] * norm_hess + tmp_mat + tmp_mat.transpose() + norm *
        //                        v_hat_deriv.transpose() * hess_hat[i] * v_hat_deriv; for (int l = 0; l < dim; l++) {
        //							(*hess)[i] += norm * deriv_hat(i, l) * v_hat_hess[l];
        //						}
        //					}
        //				}
        //            );
        for (int i = 0; i < coeff.size(); i++) {
            Eigen::MatrixXd tmp_mat = norm_deriv * deriv_hat.row(i) * v_hat_deriv;
            (*hess)[i] = coeff_hat[i] * norm_hess + tmp_mat + tmp_mat.transpose() +
                         norm * v_hat_deriv.transpose() * hess_hat[i] * v_hat_deriv;

            for (int l = 0; l < dim; l++) {
                (*hess)[i] += norm * deriv_hat(i, l) * v_hat_hess[l];
            }
        }
    }

    return coeff;
}

// Test Coeffient derivatives and hessians
void TestKruskalCoefficientDerivs(const Eigen::VectorXd& v, int ord) {
    Eigen::VectorXd rand_v = Eigen::VectorXd::Random(v.size());

    Eigen::MatrixXd deriv;
    std::vector<Eigen::MatrixXd> hess;
    Eigen::VectorXd coeff = GetKruskalCoefficients(v, ord, &deriv, &hess);

    std::cout << "\n=== order: " << ord << ", vec: " << v.transpose() << std::endl;
    std::cout << "coeff " << coeff.transpose() << std::endl;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(0.1, i);
        Eigen::MatrixXd deriv1;
        Eigen::VectorXd v1 = v + eps * rand_v;
        Eigen::VectorXd coeff1 = GetKruskalCoefficients(v1, ord, &deriv1);

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
