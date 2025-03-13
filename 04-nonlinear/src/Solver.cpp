#include "../include/Solver.h"
#include <algorithm>
#include <cmath>
#include <iostream>

/**
 * @brief Constructor for the Solver class.
 */
Solver::Solver(size_t N, double dx, double k0, double sigma, double beta,
               double delta, MPI_Comm comm)
    : N(N), dx(dx), k0(k0), sigma(sigma), beta(beta), delta(delta), comm(comm),
      builder(N, dx, k0, sigma, beta, delta), hypre(N, comm) {}

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

  for (int iter = 0; iter < max_iters; iter++) {
    // Compute dt
    double u_max = *std::max_element(u.begin(), u.end());
    dt = gamma * 2.0 /
         (4 * sigma * std::pow(u_max, 3) +
          4 * k0 * std::pow(u_max, 2) / std::pow(dx, 2));

    // Build matrix and RHS
    builder.buildMatrix(A, u);
    builder.buildRHS(F, u);

    // Assemble and solve using HYPRE
    hypre.assembleMatrix(A);
    hypre.assembleRHS(F);
    hypre.finalize();
    hypre.solveSystem(u);

    // Compute error
    double local_error = 0.0;
    for (size_t i = 0; i < N; i++)
      local_error += std::pow(u[i] - 1.0, 2);

    double global_error;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm);
    global_error = std::sqrt(global_error) / N;

    if (rank == 0)
      std::cout << "[Implicit] Iteration " << iter << " dt: " << dt
                << " Error: " << global_error << std::endl;

    if (global_error < tol)
      break;
  }
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

  // Boundary conditions
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

  for (int iter = 0; iter < max_iters; iter++) {
    // Compute the residual F
    builder.buildRHS(F, u);

    // Compute the Jacobian matrix
    computeJacobian(J, u);

    // Solve the linear system J * du = -F
    hypre.assembleMatrix(J);
    hypre.assembleRHS(F);
    hypre.finalize();
    hypre.solveSystem(du);

    // Update the solution u = u + du
    for (size_t i = 0; i < N; i++)
      u[i] += du[i];

    // Compute error
    double local_error = 0.0;
    for (size_t i = 0; i < N; i++)
      local_error += std::pow(du[i], 2);

    double global_error;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm);
    global_error = std::sqrt(global_error) / N;

    if (rank == 0)
      std::cout << "[Newton] Iteration " << iter << " Error: " << global_error
                << std::endl;

    if (global_error < tol)
      break;
  }
}
