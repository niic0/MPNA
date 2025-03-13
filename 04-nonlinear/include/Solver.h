#pragma once

#include "HypreInterface.h"
#include "MatrixBuilder.h"

#include <mpi.h>
#include <vector>

/**
 * @brief Solver class for nonlinear diffusion equation.
 */
class Solver {
public:
  Solver(size_t N, double dx, double k0, double sigma, double beta,
         double delta, MPI_Comm comm);

  void solveImplicitScheme(std::vector<double> &u, double tol, int max_iters);
  void solveNewtonRaphson(std::vector<double> &u, double tol, int max_iters);

private:
  size_t N; ///< Number of grid points
  double dx, k0, sigma, beta, delta;
  MPI_Comm comm;
  MatrixBuilder builder;
  HypreInterface hypre;

  void computeJacobian(std::vector<std::vector<double>> &J,
                       const std::vector<double> &u);
};