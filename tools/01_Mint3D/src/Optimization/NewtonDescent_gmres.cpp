
#include <fstream>
#include <iomanip>

#include <Eigen/CholmodSupport>

#include "LineSearch.h"
#include "NewtonDescent.h"
#include "CholeskyPreconditioner.h"
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
// #include <Eigen/Unsupported>
#include <Eigen/Sparse>

#include "Timer.h"
/*
namespace OptSolver {
	// Newton solver with line search


	// Eigen::GMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner> getSolverGMRES(Eigen::SparseMatrix<double>& A) {
	// 	// Create the preconditioner object
	// 	CholeskyPreconditioner preconditioner(A);

	// 	// Ensure decomposition is performed before solving
	// 	preconditioner.decompose();

	// 	// Check if decomposition was successful
	// 	if (!preconditioner.decomposed()) {
	// 		// Handle decomposition failure (e.g., use a different preconditioner or exit)
	// 		std::cerr << "Error: Cholesky decomposition failed. Aborting!" << std::endl;
	// 		throw std::runtime_error("Cholesky decomposition failed!");
	// 		// return Eigen::GMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner>(A); // Or throw an exception
	// 	}

	// 	// Create GMRES solver with the preconditioner
	// 	Eigen::GMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner> solver(A);

	// 	// Set solver options (e.g., maximum iterations, tolerance)
	// 	solver.setMaxIterations(100);
	// 	solver.setTolerance(1e-12);

	// 	return solver;
	// }

	void NewtonSolver(std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> obj_func, std::function<double(const Eigen::VectorXd&)> log_obj_func, std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> find_max_step, Eigen::VectorXd& x0, double& reg, int num_iter, double grad_tol, double xTol, double fTol, bool is_proj, bool display_info, SolveData * solve_data) {
		const int DIM = x0.rows();
		//Eigen::VectorXd randomVec = x0;
		//randomVec.setRandom();
		//x0 += 1e-6 * randomVec;
		Eigen::VectorXd grad = Eigen::VectorXd::Zero(DIM);
		Eigen::SparseMatrix<double> hessian;

		Eigen::VectorXd neg_grad, delta_x;
		double max_step_size = 1.0;
		double reg_inp = reg;
		// double reg = 1e-8;

		// is_proj = !is_proj;
		Timer<std::chrono::high_resolution_clock> total_timer;
		double total_assembling_time = 0;
		double total_solving_time = 0;
		double total_linesearch_time = 0;

		total_timer.start();
		std::ofstream optInfo;
		if (display_info) {
			std::cout << "gradient tol: " << grad_tol << ", function update tol: " << fTol << ", variable update tol: " << xTol << ", maximum iteration: " << num_iter << std::endl << std::endl;
		}
		int i = 0;

		double f = obj_func(x0, nullptr, nullptr, false);
		if (f == 0) {
			std::cout << "energy = 0, return" << std::endl;
		}

		Eigen::SparseMatrix<double> I(DIM, DIM);
		I.setIdentity();

		bool is_small_perturb_needed = false;

		double max_grad_prev = 1;
		double max_grad_curr = 1;

		double alpha_line_search = 1.0;
		double beta_kinetic_diagonal_reg = 1.0 - 1e-7;
		double kinetic_term_weight = 100;

		int toggle = 0;
		int toggle_max = 6;

		for (int i = 0; i < num_iter; i++)
		{
			if (display_info) {
				std::cout << "\niter: " << i << std::endl;
			}

			Timer<std::chrono::high_resolution_clock> local_timer;
			local_timer.start();
			double f = obj_func(x0, &grad, &hessian, is_proj);
			local_timer.stop();
			double localAssTime = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			total_assembling_time += localAssTime;

			local_timer.start();
			Eigen::SparseMatrix<double> H = hessian;
			std::cout << "num of nonzeros: " << H.nonZeros() << ", rows: " << H.rows() << ", cols: " << H.cols() << ", Sparsity: " << H.nonZeros() * 100.0 / (H.rows() * H.cols()) << "%" << std::endl;


			bool using_actual_hessian = false;
			if (is_small_perturb_needed && is_proj) { // is_small_perturb_needed && 
				// due to the numerical issue, we may need to add a small perturbation to the PSD projected hessian matrix
				H += 1e-10 * I;
				std::cout << "using proj hessian" << std::endl;
			}
			else 
			{
				H += 1e-16 * I;
				std::cout << "using un-regularized hessian" << std::endl;
				using_actual_hessian = true;
			}

			Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver(H);

			if (solver.info() != Eigen::Success)
			{
				using_actual_hessian = false;
				double b = beta_kinetic_diagonal_reg;
				double kinetic_term = (1. - b) * (1. - b) / (b * b); 
				H += (reg + kinetic_term_weight * kinetic_term) * I;
				std::cout << "using diagonally regualized hessian" << std::endl;
				Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver(H);
			}

			if ( using_actual_hessian )
			{
				if ( alpha_line_search > .9 && grad.norm() < 1e-10 )
				{
					break;
				}
			}

			//Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver(H);

			//Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(H);

			bool increase_reg = false;
			bool switched = false;



			while (solver.info() != Eigen::Success) {
				if (display_info)
				{
					if (is_proj) {
						std::cout << "some small perturb is needed to remove round-off error, current reg = " << reg << std::endl;
					}

					else {
						std::cout << "Matrix is not positive definite, current reg = " << reg <<  " | current kinetic term " << std::pow(beta_kinetic_diagonal_reg,2) <<  std::endl;
					}

				}

				if (is_proj) {
					is_small_perturb_needed = true;
				}

				increase_reg = true;


				// H = hessian + reg * I;

				double b = beta_kinetic_diagonal_reg;
				double kinetic_term = (1. - b) * (1. - b) / (b * b); 
				H = hessian + (reg + kinetic_term_weight * kinetic_term) * I;

				solver.compute(H);
				reg = std::max(2 * reg, 1e-16);
				beta_kinetic_diagonal_reg = std::max( beta_kinetic_diagonal_reg * 0.5, 1e-8 );


				if (reg > 1e12) {
					std::cout << "reg is too large, switched to using SPD hessian instead." << std::endl;
					reg = 1e-8;
					is_proj = true;
					f = obj_func(x0, &grad, &hessian, is_proj);
					H = hessian + 1e-8 * I;
					solver.compute(H);
				}
			}

			neg_grad = -grad;


			CholeskyPreconditioner P(H);

			Eigen::GMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner> gmres;
			// Eigen::DGMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner> gmres;

			
			gmres.compute(H);
			gmres.setTolerance(1e-12);
			gmres.setMaxIterations(1000);

			// gmres.setMaxIterations(500);
			// gmres.set_restart(30);
			// gmres.setEigenv(10);  // Number of deflated eigenvalues

			// gmres.set_restart(100);
			// gmres.setEigenv(30);


			// delta_x = gmres.solve(neg_grad); 

			delta_x = gmres.solveWithGuess(neg_grad, neg_grad); 


			// Eigen::GMRES<Eigen::SparseMatrix<double>, CholeskyPreconditioner> gmres;
			// gmres.compute(H);
			// gmres.setPreconditioner(P);
			// delta_x = gmres.solve(neg_grad);

			// Output results
			// std::cout << "The solution is:\n" << delta_x << std::endl;
			std::cout << "Number of iterations: " << gmres.iterations() << std::endl;
			std::cout << "Estimated error: " << gmres.error() << std::endl;

			// auto solver_gmres = getSolverGMRES(H);
			// delta_x = solver.solveWithGuess(neg_grad, neg_grad);


			// delta_x = solver.solve(neg_grad);

			local_timer.stop();
			double local_solving_time = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			total_solving_time += local_solving_time;


			max_step_size = find_max_step(x0, delta_x);

			local_timer.start();
			double rate = BacktrackingArmijo(x0, grad, delta_x, obj_func, max_step_size);
			local_timer.stop();
			double local_linesearch_time = local_timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			total_linesearch_time += local_linesearch_time;
			alpha_line_search = rate;


			// If line search progress is small, toggle next step projection 
			if ( rate < 1e-3 )
			{
				
				toggle = (toggle + 1) % toggle_max;

				if ( toggle == 0 )
				{
					is_proj = !is_proj;
					switched = true;
				}
			}

			// Update kinetic term and exponential term.  
			if (!is_proj) {
				// Update the control parameter of the kinetic term, scales like 1 / (beta)^2
				if (rate < .3)
				{
					beta_kinetic_diagonal_reg = std::max( beta_kinetic_diagonal_reg * 0.5, 1e-8 );
				}
				else if (rate > .9)
				{
					beta_kinetic_diagonal_reg = std::min( beta_kinetic_diagonal_reg * 2, 1 - 1e-7 );
				}

				reg *= 0.8;
				reg = std::max(reg, 1e-16);
				reg = std::min(reg, 1e16);

				// if ( rate < .001 )
				// {
				// 	reg *= 10;
				// }

				// if ( reg > 1e12)
				// {
				// 	reg = 1e-12;
				// }
			
			}
			else {
				// Exit condition for projection state.  
				if ( (((rate < .4 ) && grad.norm() < .005) || ( rate > .999 &&  grad.norm() < .001 )) && !switched )// delta_x.norm() < 1e-8) )
				{
					is_proj = false;
					reg = 1e-12;
					std::cout << "Switch back to diag hessian " << reg << std::endl;
				}
				// reg = 1e-8;
			}

			x0 = x0 + rate * delta_x;

			double fnew = obj_func(x0, &grad, nullptr, is_proj);
			double residual = (H * delta_x + grad).norm();
			log_obj_func(x0);
			if (display_info) {
				std::cout << "line search rate : " << rate << ", actual hessian : " << !is_proj << ", reg = " << reg << std::endl;
				std::cout << "f_old: " << f << ", f_new: " << fnew << ", grad norm: " << grad.norm() << ", delta x: " << rate * delta_x.norm() << ", delta_f: " << f - fnew  << ", residual: " << residual << std::endl;
				std::cout << "timing info (in total seconds): " << std::endl;
				std::cout << "assembling took: " << total_assembling_time << ", LLT solver took: " << total_solving_time << ", line search took: " << total_linesearch_time << std::endl;
			}

			double b = beta_kinetic_diagonal_reg;
			double kinetic_term = (1. - b) * (1. - b) / (b * b); 

			solve_data->total_time.push_back(total_timer.elapsed<std::chrono::milliseconds>() * 1e-3);
			solve_data->assembly_time.push_back(total_assembling_time);
			solve_data->solve_time.push_back(total_solving_time);
			solve_data->line_search_time.push_back(total_linesearch_time);
			solve_data->identity_weight.push_back(reg + kinetic_term_weight * kinetic_term);
			solve_data->total_energy.push_back(fnew);
			solve_data->energy_diff.push_back(f - fnew);
			solve_data->solve_residual.push_back( residual );
			solve_data->rhs_norm.push_back( grad.norm() );


			// // switch to the actual hessian when close to convergence
			// if ((f - fnew) / f < 1e-5 || delta_x.norm() < 1e-5 || grad.norm() < 1e-4) {
			// 	is_proj = false;
			// }



			is_proj = is_proj || increase_reg;


			// // Termination conditions
			// if (rate < 1e-8) {
			// 	std::cout << "terminate with small line search rate (<1e-8): L2-norm = " << grad.norm() << std::endl;
			// 	break;
			// }

			if (grad.norm() < grad_tol && rate > .999) {
				std::cout << "terminate with gradient L2-norm = " << grad.norm() << std::endl;
				break;
			}

			// if (rate * delta_x.norm() < xTol && rate > .999) {
			// 	std::cout << "terminate with small variable change, gradient L2-norm = " << grad.norm() << std::endl;
			// 	break;
			// }

			// if (f - fnew < fTol) {
			// 	std::cout << "terminate with small energy change, gradient L2-norm = " << grad.norm() << std::endl;
			// 	break;
			// }
		}




		if (i >= num_iter) {
			std::cout << "terminate with reaching the maximum iteration, with gradient L2-norm = " << grad.norm() << std::endl;
		}

		f = obj_func(x0, &grad, nullptr, false);
		std::cout << "end up with energy: " << f << ", gradient: " << grad.norm() << std::endl;

		total_timer.stop();
		if (display_info) {
			std::cout << "total time cost (s): " << total_timer.elapsed<std::chrono::milliseconds>() * 1e-3 << ", within that, assembling took: " << total_assembling_time << ", LLT solver took: " << total_solving_time << ", line search took: " << total_linesearch_time << std::endl;
		}



	}

	// Test the function gradient and hessian
	void TestFuncGradHessian(std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> obj_Func, const Eigen::VectorXd& x0)
	{
		Eigen::VectorXd dir = x0;
		dir(0) = 0;
		dir.setRandom();

		Eigen::VectorXd grad;
		Eigen::SparseMatrix<double> H;

		double f = obj_Func(x0, &grad, &H, false);
		std::cout << "f: " << f << std::endl;
		if (f == 0)
			return;

		for (int i = 3; i < 10; i++)
		{
			double eps = std::pow(0.1, i);
			Eigen::VectorXd x = x0 + eps * dir;
			Eigen::VectorXd grad1;
			double f1 = obj_Func(x, &grad1, nullptr, false);

			std::cout << "\neps: " << eps << std::endl;
			std::cout << "energy-gradient: " << (f1 - f) / eps - grad.dot(dir) << std::endl;
			std::cout << "gradient-hessian: " << ((grad1 - grad) / eps - H * dir).norm() << std::endl;
		}
	}
}


*/