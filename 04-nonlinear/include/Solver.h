#ifndef SOLVER_H
#define SOLVER_H

#include "HypreInterface.h"
#include "MatrixBuilder.h"
#include <chrono>
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
  void printSummary() const;

private:
  size_t N;
  double dx, k0, sigma, beta, delta;
  MPI_Comm comm;
  MatrixBuilder builder;
  HypreInterface hypre;

  int implicit_iterations;
  int newton_iterations;
  double implicit_time;
  double newton_time;
  double implicit_error;
  double newton_error;

  void computeJacobian(std::vector<std::vector<double>> &J,
                       const std::vector<double> &u);
};

#endif // SOLVER_H
