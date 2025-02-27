#pragma once

#include <cstddef>
#include "CSRMatrix.h"

/**
 * @brief Structure to store solver performance results.
 */
struct SolverResult {
  std::string solver_name;
  int iterations;
  double residual;
  double time_ms;
};

size_t toMatId(size_t i, size_t j, size_t mesh_size);
double computeResidual(CSRMatrix &A, const std::vector<double> &x,
                       const std::vector<double> &b);
void printResultsTable(const std::vector<SolverResult> &results);