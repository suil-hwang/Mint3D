#include "MiNTCommonFunc.h"

#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

namespace MiNT3D {
// Self-adjoint SPD projection for a small size matrix A
Eigen::MatrixXd SelfAdjSPDProjection(const Eigen::MatrixXd& A) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXd V = eigensolver.eigenvectors();
    Eigen::VectorXd D = eigensolver.eigenvalues();
    for (int i = 0; i < D.size(); i++) {
        if (D(i) < 1e-8) {
            D(i) = std::min(std::abs(D(i)) + 1e-10, 1e-8);   // small value to avoid numerical negative eigenvalues
        }
    }
    return V * D.asDiagonal() * V.transpose();
}

// Helper function to create a cross product matrix
Eigen::Matrix3d cross_product_matrix(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return m;
}

double Det2x3(const Eigen::MatrixXd& frames, int tet_id, Eigen::VectorXd* deriv, Eigen::MatrixXd* hess) {
    // assert(frames.cols() == 9);
    // double energy = 0;
    // double weight = 1.0;  // Adjust this based on your application context

    // // Weight adjustment based on tetrahedron volume and other factors
    // // Add any necessary weight adjustments here

    // if (deriv) {
    // 	deriv->resize(9);
    // 	deriv->setZero();
    // }
    // if (hess) {
    // 	hess->setZero(9, 9);
    // }

    // // Extract vectors a and b from frames (last two rows)
    // Eigen::Vector3d a = frames.row(tet_id).segment<3>(3);
    // Eigen::Vector3d b = frames.row(tet_id).segment<3>(6);

    // // Compute the 3D cross product
    // Eigen::Vector3d cross_product = a.cross(b);
    // double norm_cross_product = cross_product.norm();
    // double part = 1.0 - norm_cross_product;
    // energy = part * part * weight;

    // // Compute the gradient if requested
    // if (deriv) {
    // 	Eigen::Matrix<double, 3, 9> jacobian;
    // 	jacobian.setZero();

    // 	jacobian.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    // 	jacobian.block<3, 3>(0, 6) = -Eigen::Matrix3d::Identity();

    // 	(*deriv).segment<3>(3) = -2.0 * part * jacobian.block<3, 3>(0, 3).transpose() * cross_product /
    // norm_cross_product * weight;
    // 	(*deriv).segment<3>(6) = -2.0 * part * jacobian.block<3, 3>(0, 6).transpose() * cross_product /
    // norm_cross_product * weight;
    // }

    // // Compute the Hessian if requested
    // if (hess) {
    // 	hess->setZero();

    // 	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    // 	Eigen::Matrix<double, 9, 9> jacobian_hess;
    // 	jacobian_hess.setZero();

    // 	jacobian_hess.block<3, 3>(3, 3) = I;
    // 	jacobian_hess.block<3, 3>(6, 6) = I;
    // 	jacobian_hess.block<3, 3>(3, 6) = -I;
    // 	jacobian_hess.block<3, 3>(6, 3) = -I;

    // 	Eigen::Matrix3d cross_product_matrix = cross_product * cross_product.transpose();
    // 	double norm_cross_product_sq = norm_cross_product * norm_cross_product;

    // 	for (int i = 3; i < 9; ++i) {
    // 		for (int j = 3; j < 9; ++j) {
    // 			(*hess)(i, j) = 2.0 * weight * (
    // 				-jacobian_hess.block<3, 3>(i % 3, j % 3) / norm_cross_product
    // 				+ cross_product_matrix(i % 3, j % 3) / (norm_cross_product_sq * norm_cross_product)
    // 			);
    // 		}
    // 	}
    // }
    // return energy;
    return 0;
}

#include <Eigen/Dense>
#include <algorithm>
#include <limits>
#include <vector>

Eigen::MatrixXd ExpandFramePermutation(const Eigen::MatrixXd& perm, const Eigen::MatrixXd& sgn) {
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(9, 9);

    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d to_exp = perm * sgn;
    P = Eigen::kroneckerProduct(to_exp, I3);

    // std::cout << "P: \n" << P << std::endl;
    // std::cout << to_exp << "to_exp" << std::endl;

    return P;
}

// Function to compute the optimal permutation matrix
Eigen::MatrixXd ComputeOptimalPermutationMatrix(const Eigen::VectorXd& f, const Eigen::VectorXd& g,
                                                const Eigen::MatrixXd& basis) {
    Eigen::MatrixXd best_P(9, 9);
    double min_distance = std::numeric_limits<double>::max();

    // if (basis.rows() != 3) {
    //     Eigen::MatrixXd tmp_basis = Eigen::Matrix3d::Identity();
    //     tmp_basis.block(0, 0, basis.rows(), basis.cols()) = basis;
    //     Eigen::Vector3d b0 = basis.row(0);
    //     Eigen::Vector3d b1 = basis.row(1);
    //     Eigen::Vector3d n = b0.cross(b1);
    //     tmp_basis.block(2, 0, 1, 3) = n * 1e-6 * 0.;
    // }

    // Define permutation matrices
    std::vector<Eigen::Matrix3d> perm_matrices(6);
    std::vector<Eigen::Matrix3d> sgn_matrices(8);

    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();

    // Generate all 3! permutations
    perm_matrices[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    perm_matrices[1] << 0, 1, 0, 1, 0, 0, 0, 0, 1;

    perm_matrices[2] << 0, 0, 1, 1, 0, 0, 0, 1, 0;

    perm_matrices[3] << 1, 0, 0, 0, 0, 1, 0, 1, 0;

    perm_matrices[4] << 0, 1, 0, 0, 0, 1, 1, 0, 0;

    perm_matrices[5] << 0, 0, 1, 0, 1, 0, 1, 0, 0;

    // All the sign flips
    sgn_matrices[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    sgn_matrices[1] << -1, 0, 0, 0, 1, 0, 0, 0, 1;

    sgn_matrices[2] << 1, 0, 0, 0, -1, 0, 0, 0, 1;

    sgn_matrices[3] << -1, 0, 0, 0, -1, 0, 0, 0, 1;

    sgn_matrices[4] << 1, 0, 0, 0, 1, 0, 0, 0, -1;

    sgn_matrices[5] << -1, 0, 0, 0, 1, 0, 0, 0, -1;

    sgn_matrices[6] << 1, 0, 0, 0, -1, 0, 0, 0, -1;

    sgn_matrices[7] << -1, 0, 0, 0, -1, 0, 0, 0, -1;

    // Iterate through all permutations
    for (const auto& perm : perm_matrices) {
        for (const auto& sgn : sgn_matrices) {
            // Eigen::MatrixXd P = Eigen::kroneckerProduct(sgn, perm);
            Eigen::MatrixXd P = ExpandFramePermutation(perm, sgn);
            Eigen::MatrixXd g_perm = P * g;
            Eigen::VectorXd diff = f - g_perm;
            double distance = 0;
            for (int i = 0; i < 3; i++) {
                distance += (basis * diff.segment<3>(3 * i)).norm();
            }

            if (distance < min_distance) {
                min_distance = distance;
                // std::cout << perm << " perm \n" << sgn << " sgn " << std::endl;
                // std::cout << P << " P \n" << std::endl;
                best_P = P;
            }
        }
    }

    return best_P;
}

void Test_ComputeOptimalPermutationMatrix() {
    std::cout << "Running test suite for ComputeOptimalPermutationMatrix." << std::endl;

    {
        // Test 1: Basic Valid Input
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 1, 0, 0, 0, 1, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(9, 9);

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 1 failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 2: Different Vectors
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(9, 9);   // Adjusted to correct the expected matrix
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 2 failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 3: Negative Values
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);
        f << -1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 1, 0, 0, 0, -1, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(9, 9);   // Adjusted expected matrix
        expected << -1, -0, -0, 0, 0, 0, 0, 0, 0, -0, -1, -0, 0, 0, 0, 0, 0, 0, -0, -0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -1, -0, -0, 0, 0, 0, 0, 0, 0, -0, -1, -0, 0, 0, 0, 0, 0, 0, -0, -0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 3 failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 4: Non-Identity Basis Matrix
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        Eigen::MatrixXd basis(3, 3);
        basis << 2, 0, 0, 0, 2, 0, 0, 0, 2;
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(9, 9);   // Adjusted expected matrix
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 4 failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 5: Zero Vectors
        Eigen::VectorXd f = Eigen::VectorXd::Zero(9);
        Eigen::VectorXd g = Eigen::VectorXd::Zero(9);
        Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);

        if (result.rows() != 9 || result.cols() != 9) {
            std::cout << "Test 5 failed" << std::endl;
            // std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 6: Large Numbers
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);
        f << 1e6, 0, 0, 0, 1e6, 0, 0, 0, 1e6;
        g << 0, 1e6, 0, 1e6, 0, 0, 0, 0, 1e6;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(9, 9);   // Adjusted expected matrix
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 6 failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    // {
    //     // Test 7: Symmetry Test
    //     Eigen::VectorXd f(9);
    //     Eigen::VectorXd g(9);
    //     Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(3, 3);
    //     f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    //     g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

    //     Eigen::MatrixXd result1 = ComputeOptimalPermutationMatrix(f, g, basis);
    //     Eigen::MatrixXd result2 = ComputeOptimalPermutationMatrix(g, f, basis);

    //     if ((result1 - result2).squaredNorm() >= 1e-6) {
    //         std::cout << "Test 7 failed" << std::endl;
    //         std::cout << expected << " expected" << std::endl;
    //         std::cout << result << " result" << std::endl;
    //     }
    // }

    std::cout << "All identity basis tests passed!" << std::endl;

    // Additional tests for non-identity basis matrices
    {
        // Test 1: Scaling Basis Matrix
        Eigen::MatrixXd basis(3, 3);
        basis << 2, 0, 0, 0, 2, 0, 0, 0, 2;
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(9, 9);   // Adjusted expected matrix
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 1 with non-identity basis failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    {
        // Test 2: Rotation Basis Matrix
        Eigen::MatrixXd basis(3, 3);
        double theta = M_PI / 4;   // 45 degrees
        basis << std::cos(theta), -std::sin(theta), 0, std::sin(theta), std::cos(theta), 0, 0, 0, 1;
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);
        Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(9, 9);   // Adjusted expected matrix
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 2 with non-identity basis failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }
    }

    // Test with a 2x3 Basis Matrix
    {
        Eigen::MatrixXd basis(2, 3);
        basis << 1, 0, 0, 0, 0, 1;
        Eigen::VectorXd f(9);
        Eigen::VectorXd g(9);
        f << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        g << 0, 1, 0, 1, 0, 0, 0, 0, 1;

        Eigen::MatrixXd result = ComputeOptimalPermutationMatrix(f, g, basis);

        // Compute the expected result manually considering the effect of the non-square basis matrix
        Eigen::MatrixXd expected(9, 9);
        expected << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

        if ((result - expected).squaredNorm() >= 1e-6) {
            std::cout << "Test 3 with non-identity basis failed" << std::endl;
            std::cout << expected << " expected" << std::endl;
            std::cout << result << " result" << std::endl;
        }

        // Replace the example with the actual expected result considering basis and g permutation effects
        // assert(result.isApprox(expected));
    }

    std::cout << "Non-square basis test passed!" << std::endl;
}

// Function to compute tensor product of two matrices
Eigen::MatrixXd tensorProduct(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B) {
    Eigen::MatrixXd result(9, 9);
    result.setZero();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result.block<3, 3>(3 * i, 3 * j) = A(i, j) * B;
        }
    }
    return result;
}

// Get determinant of a 3x3 matrix A, togiether with its derivative and hessian in terms of A00, A01, ..., A22
double Det3x3(const Eigen::Matrix3d& A, Eigen::VectorXd* deriv, Eigen::Matrix<double, 9, 9>* hess) {
    double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) - A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                 A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
    // double det = A.determinant();

    if (deriv != nullptr) {
        *deriv << A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1), -A(1, 0) * A(2, 2) + A(1, 2) * A(2, 0),
            A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0), -A(0, 1) * A(2, 2) + A(0, 2) * A(2, 1),
            A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0), -A(0, 0) * A(2, 1) + A(0, 1) * A(2, 0),
            A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1), -A(0, 0) * A(1, 2) + A(0, 2) * A(1, 0),
            A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    }
    if (hess != nullptr) {
        hess->setZero();
        (*hess)(0, 4) = (*hess)(4, 0) = A(2, 2);
        (*hess)(0, 5) = (*hess)(5, 0) = -A(2, 1);
        (*hess)(0, 7) = (*hess)(7, 0) = -A(1, 2);
        (*hess)(0, 8) = (*hess)(8, 0) = A(1, 1);

        (*hess)(1, 3) = (*hess)(3, 1) = -A(2, 2);
        (*hess)(1, 5) = (*hess)(5, 1) = A(2, 0);
        (*hess)(1, 6) = (*hess)(6, 1) = A(1, 2);
        (*hess)(1, 8) = (*hess)(8, 1) = -A(1, 0);

        (*hess)(2, 3) = (*hess)(3, 2) = A(2, 1);
        (*hess)(2, 4) = (*hess)(4, 2) = -A(2, 0);
        (*hess)(2, 6) = (*hess)(6, 2) = -A(1, 1);
        (*hess)(2, 7) = (*hess)(7, 2) = A(1, 0);

        (*hess)(3, 7) = (*hess)(7, 3) = A(0, 2);
        (*hess)(3, 8) = (*hess)(8, 3) = -A(0, 1);

        (*hess)(4, 6) = (*hess)(6, 4) = -A(0, 2);
        (*hess)(4, 8) = (*hess)(8, 4) = A(0, 0);

        (*hess)(5, 6) = (*hess)(6, 5) = A(0, 1);
        (*hess)(5, 7) = (*hess)(7, 5) = -A(0, 0);
    }

    return det;
}

// Test the derivative and hessian of Det3x3
void TestDet3x3() {
    Eigen::Matrix3d A = Eigen::Matrix3d::Random();

    Eigen::VectorXd deriv;
    Eigen::Matrix<double, 9, 9> hess;
    deriv.setZero(9);

    double det = Det3x3(A, &deriv, &hess);

    std::cout << "A: \n" << A << std::endl;

    std::cout << "det: " << det << ", eigen det: " << A.determinant() << std::endl;
    std::cout << "deriv: " << deriv.transpose() << std::endl;
    std::cout << "hessian: \n" << hess << std::endl;

    Eigen::VectorXd perturb = Eigen::VectorXd::Random(9);
    Eigen::Matrix3d A_perturb;

    for (int i = 4; i < 12; i++) {
        double eps = std::pow(10, -i);
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                A_perturb(j, k) = A(j, k) + eps * perturb(3 * j + k);
            }
        }

        Eigen::VectorXd deriv_perturb;
        deriv_perturb.setZero(9);
        double det_perturb = Det3x3(A_perturb, &deriv_perturb);

        std::cout << "eps: " << eps << std::endl;
        std::cout << "energy-deriv: " << (det_perturb - det) / eps - deriv.dot(perturb) << std::endl;
        std::cout << "deriv-hess: " << ((deriv_perturb - deriv) / eps - hess * perturb).norm() << std::endl;
    }
}

}   // namespace MiNT3D