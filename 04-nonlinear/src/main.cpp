#include "../include/HypreInterface.h"
#include "../include/MatrixBuilder.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

/**
 * @brief Computes the thermal conductivity k(u).
 * @param u The dimensionless temperature variable u = T / T0.
 * @param k0 The reference thermal conductivity.
 * @return double The computed value of κ(u).
 */
double k(double u, double k0) { return k0 * u * u; }

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Set up parameters
  size_t N = 10000; // Number of grid points
  double L = 1.0;
  double dx = L / (N - 1);
  double k0 = 0.01, sigma = 0.1, beta = 1.0, delta = 0.2;
  double gamma = 0.1, dt;
  double tol = 1e-5;

  // Initialize the problem
  MatrixBuilder builder(N, dx, k0, sigma, beta, delta);
  HypreInterface hypre(N, MPI_COMM_WORLD);

  // Create matrices and vectors
  std::vector<std::vector<double>> A;
  std::vector<double> F;
  std::vector<double> u(N, 1.0); // Initial temperature field
  u[0] = 1.1;                    // Ensures u_max > 1

  // Compute the discretization
  builder.buildMatrix(A, u);
  builder.buildRHS(F, u);

  int max_iters = 1000;
  for (int iter = 0; iter < max_iters; iter++) {
    // Synchronize u across all processes before computing dt
    std::vector<double> global_u(N);
    MPI_Allreduce(u.data(), global_u.data(), N, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    // Compute u_max correctly
    double local_u_max = *std::max_element(global_u.begin(), global_u.end());
    double u_max;
    MPI_Allreduce(&local_u_max, &u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double denominator =
        4 * sigma * pow(u_max, 3) + (4 * k(u_max, k0)) / pow(dx, 2);

    if (denominator < 1e-10 || std::isnan(denominator) ||
        std::isinf(denominator)) {
      denominator = 1e-10; // Prevent division by zero
    }

    dt = gamma * 2.0 / denominator;

    if (dt < 1e-10 || std::isnan(dt) || std::isinf(dt)) {
      if (rank == 0) {
        std::cerr << "⚠️  ERROR: dt is too small or NaN/inf. dt = " << dt
                  << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Build matrix and RHS
    builder.buildMatrix(A, u);
    builder.buildRHS(F, u);

    // Assemble matrix and RHS in HYPRE
    hypre.assembleMatrix(A);
    hypre.assembleRHS(F);
    hypre.finalize();

    // Solve system
    hypre.solveSystem(u);

    // Synchronize u again to ensure all processes have updated values
    MPI_Allreduce(u.data(), global_u.data(), N, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    u = global_u; // Copy back to local u

    // Compute error
    double local_error = 0.0;
    for (size_t i = 0; i < N; i++)
      local_error += std::pow(u[i] - 1.0, 2);

    double global_error;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    global_error = std::sqrt(global_error) / N;

    if (rank == 0)
      std::cout << "Iteration " << iter << " u_max: " << u_max << " dt: " << dt
                << " Error: " << global_error << std::endl;

    if (global_error < tol)
      break;
  }

  MPI_Finalize();
  return 0;
}
