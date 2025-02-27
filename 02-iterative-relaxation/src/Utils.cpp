#include <vector>
#include <cmath>
#include <iomanip>

#include "../include/CSRMatrix.h"
#include "../include/Utils.h"

double computeResidual(CSRMatrix &A, const std::vector<double> &x,
                       const std::vector<double> &b) {
  double norm = 0.0;
  std::vector<double> Ax(A.getRows(), 0.0);
  A.matrixVectorMultiply(Ax, x); // Compute Ax product

  for (size_t i = 0; i < A.getRows(); i++) {
    double r = b[i] - Ax[i];
    norm += r * r;
  }

  return std::sqrt(norm);
}

/**
 * @brief Prints a formatted table of solver results.
 *
 * @param results A vector containing solver performance data.
 */
void printResultsTable(const std::vector<SolverResult>& results) {
    std::cout << std::left << std::setw(20) << "Solver" 
              << std::setw(20) << "Iterations" 
              << std::setw(20) << "Residual" 
              << std::setw(20) << "Time (ms)" << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (const auto& res : results) {
        std::cout << std::left << std::setw(20) << res.solver_name
                  << std::setw(20) << (res.iterations == -1 ? "Not converged" : std::to_string(res.iterations))
                  << std::setw(20) << res.residual
                  << std::setw(20) << res.time_ms << "\n";
    }
}