#include "../include/Solver.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

/**
 * @brief Constructor for the Solver class.
 */
Solver::Solver(size_t N, double dx, double k0, double sigma, double beta,
               double delta, MPI_Comm comm)
    : N(N), dx(dx), k0(k0), sigma(sigma), beta(beta), delta(delta), comm(comm),
      builder(N, dx, k0, sigma, beta, delta), hypre(N, comm),
      implicit_error(0), newton_error(0), implicit_iterations(0), newton_iterations(0), implicit_time(0.0),
      newton_time(0.0) {}

/**
 * @brief Solves the diffusion equation using the implicit scheme.
 */
void Solver::solveImplicitScheme(std::vector<double> &u, double tol,
                                 int max_iters) {
  std::vector<std::vector<double>> A;
  std::vector<double> F;
  double gamma = 0.1, dt;
  int rank;
  MPI_Comm_rank(comm, &rank);

  auto start = std::chrono::high_resolution_clock::now();

  for (int iter = 0; iter < max_iters; iter++) {
    double u_max = *std::max_element(u.begin(), u.end());
    dt = gamma * 2.0 /
         (4 * sigma * std::pow(u_max, 3) +
          4 * k0 * std::pow(u_max, 2) / std::pow(dx, 2));

    builder.buildMatrix(A, u);
    builder.buildRHS(F, u);

    hypre.assembleMatrix(A);
    hypre.assembleRHS(F);
    hypre.finalize();
    hypre.solveSystem(u);

    double local_error = 0.0;
    for (size_t i = 0; i < N; i++)
      local_error += std::pow(u[i] - 1.0, 2);

    double global_error;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm);
    global_error = std::sqrt(global_error) / N;

    implicit_iterations = iter;
    implicit_error = global_error;

    if (global_error < tol) {
      break;
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  implicit_time = std::chrono::duration<double>(end - start).count();

  if (rank == 0)
    std::cout << "Implicit DONE..." << std::endl;
}

/**
 * @brief Computes the Jacobian matrix for Newton-Raphson method.
 */
void Solver::computeJacobian(std::vector<std::vector<double>> &J,
                             const std::vector<double> &u) {
  J.assign(N, std::vector<double>(N, 0.0));
  double dx2 = dx * dx;

  for (size_t i = 1; i < N - 1; i++) {
    double kappa_ip = (builder.kappa(u[i + 1]) + builder.kappa(u[i])) / 2.0;
    double kappa_im = (builder.kappa(u[i - 1]) + builder.kappa(u[i])) / 2.0;
    double dFdUi = 4.0 * sigma * std::pow(u[i], 3);

    J[i][i - 1] = -kappa_im / dx2;
    J[i][i] = (kappa_im + kappa_ip) / dx2 + dFdUi;
    J[i][i + 1] = -kappa_ip / dx2;
  }

  J[0][0] = 1.0;
  J[N - 1][N - 1] = 1.0;
}

/**
 * @brief Solves the nonlinear diffusion equation using Newton-Raphson method.
 */
void Solver::solveNewtonRaphson(std::vector<double> &u, double tol,
                                int max_iters) {
  std::vector<std::vector<double>> J;
  std::vector<double> F, du(N);
  int rank;
  MPI_Comm_rank(comm, &rank);

  auto start = std::chrono::high_resolution_clock::now();

  for (int iter = 0; iter < max_iters; iter++) {
    builder.buildRHS(F, u);
    computeJacobian(J, u);

    hypre.assembleMatrix(J);
    hypre.assembleRHS(F);
    hypre.finalize();
    hypre.solveSystem(du);

    for (size_t i = 0; i < N; i++)
      u[i] += du[i];

    double local_error = 0.0;
    for (size_t i = 0; i < N; i++)
      local_error += std::pow(du[i], 2);

    double global_error;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm);
    global_error = std::sqrt(global_error) / N;

    newton_iterations = iter;
    newton_error = global_error;

    if (global_error < tol) {
      break;
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  newton_time = std::chrono::duration<double>(end - start).count();

  if (rank == 0)
    std::cout << "Newton DONE..." << std::endl;
}

/**
 * @brief Prints a summary of the results.
 */
void Solver::printSummary() const {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    std::cout << "\n===== SOLVER SUMMARY =====\n";
    std::cout << " Implicit Scheme:\n";
    std::cout << "   - Iterations: " << implicit_iterations << "\n";
    std::cout << "   - Time: " << implicit_time << " sec\n";
    std::cout << "   - Error: " << implicit_error << "\n";
    std::cout << " Newton-Raphson:\n";
    std::cout << "   - Iterations: " << newton_iterations << "\n";
    std::cout << "   - Time: " << newton_time << " sec\n";
    std::cout << "   - Time: " << newton_error << "\n";
    std::cout << "==========================\n";
  }
}
