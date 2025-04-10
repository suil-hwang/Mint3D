#pragma once

#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

class CholeskyPreconditioner {
   public:
    using MatrixType = Eigen::SparseMatrix<double>;

    CholeskyPreconditioner() : m_info(Eigen::Success) {}   // Default constructor

    explicit CholeskyPreconditioner(const MatrixType& matrix) : m_info(Eigen::Success) { compute(matrix); }

    void compute(const MatrixType& matrix) {
        cholmod.compute(matrix);
        m_info = cholmod.info();
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const {
        return cholmod.solve(rhs);   // This needs to return a vector directly
    }

    Eigen::ComputationInfo info() const { return m_info; }

   private:
    Eigen::CholmodSupernodalLLT<MatrixType> cholmod;
    Eigen::ComputationInfo m_info;
};

// // Preconditioner class with Cholesky decomposition storage and error handling
// class CholeskyPreconditioner : public Eigen::SparseMatrix<double>::Base {
// public:
// //   EIGEN_MAKE_UFC_MFC_ mambo; // weird magic, name can be anything

//   CholeskyPreconditioner(const Eigen::SparseMatrix<double>& A) : A_(A), decomposed_(false) {}

//   // Function to perform Cholesky decomposition and store L factor
//   void decompose() {
//     // Attempt Cholesky decomposition
//     L_ = A.llt().matrixL();

//     // Check for successful decomposition (positive definite check)
//     if (L_.info() != Eigen::Success) {
//       decomposed_ = false;
//       // Handle decomposition failure (e.g., throw an exception or provide a default preconditioner)
//       std::cerr << "Error: Cholesky decomposition failed!" << std::endl;
//     } else {
//       decomposed_ = true;
//     }
//   }

//   virtual void solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) const override {
//     if (!decomposed_) {
//       // Handle the case where decomposition failed
//       return;
//     }

//     // Precondition using forward substitution with L
//     x = L_.triangularSolve<Eigen::Lower triangular, Eigen::NoDiagonalPreconditioner>(b);
//   }

// private:
//   const Eigen::SparseMatrix<double>& A_;
//   Eigen::SparseMatrix<double> L_; // Stores the L factor from Cholesky decomposition (if successful)
//   bool decomposed_;
// };
