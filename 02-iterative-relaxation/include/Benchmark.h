#pragma once

#include <chrono>
#include <functional>

#include "Solver.h"
#include "Utils.h"

SolverResult
runSolver(const std::string &solver_name,
          std::function<int(CSRMatrix &, const std::vector<double> &,
                            std::vector<double> &, size_t, double)>
              solver_function,
          CSRMatrix &A, const std::vector<double> &b,
          std::vector<double> &x, size_t max_iter, double tolerance);
