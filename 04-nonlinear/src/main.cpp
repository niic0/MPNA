#include "../include/Solver.h"

#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set up parameters
  size_t N = 10000; // Number of grid points
  double L = 1.0;
  double dx = L / (N - 1);
  double k0 = 0.01, sigma = 0.1, beta = 1.0, delta = 0.2;
  double tol = 1e-5;
  int max_iters = 10;

  // Initialize solver
  Solver solver(N, dx, k0, sigma, beta, delta, MPI_COMM_WORLD);

  // Initial condition
  std::vector<double> u(N, 1.0);
  u[0] = 1.1; // Ensures u_max > 1
  solver.solveImplicitScheme(u, tol, max_iters);

  // Reset u for Newton-Raphson
  std::fill(u.begin(), u.end(), 1.0);
  u[0] = 1.1;
  solver.solveNewtonRaphson(u, tol, max_iters);

  solver.printSummary();
  MPI_Finalize();
  return 0;
}
