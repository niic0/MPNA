#include <chrono>
#include <functional>

#include "../include/Utils.h"

/**
 * @brief Runs a solver and records its execution time and residual.
 *
 * @param solver_name The name of the solver.
 * @param solver_function The solver function to execute.
 * @param A The sparse matrix in CSR format.
 * @param b The right-hand side vector.
 * @param x The solution vector (modified in-place).
 * @param max_iter Maximum number of iterations.
 * @param tolerance Convergence tolerance.
 * @return SolverResult containing solver stats.
 */
SolverResult
runSolver(const std::string &solver_name,
          std::function<int(CSRMatrix &, const std::vector<double> &,
                            std::vector<double> &, size_t, double)>
              solver_function,
          CSRMatrix &A, const std::vector<double> &b,
          std::vector<double> &x, size_t max_iter, double tolerance) {
  std::fill(x.begin(), x.end(), 0.0); // Reset solution

  auto start = std::chrono::high_resolution_clock::now();
  int iterations = solver_function(A, b, x, max_iter, tolerance);
  auto end = std::chrono::high_resolution_clock::now();

  double residual = computeResidual(A, x, b);
  double time_ms =
      std::chrono::duration<double, std::milli>(end - start).count();

  return {solver_name, iterations, residual, time_ms};
}
