#pragma once

#include "CSRMatrix.h"

int jacobi(CSRMatrix &A, const std::vector<double> &b, std::vector<double> &x,
           size_t max_iters, double tol);

int gaussSeidel(CSRMatrix &A, const std::vector<double> &b,
                std::vector<double> &x, size_t max_iters, double tol);

int conjugateGradient(CSRMatrix &A, const std::vector<double> &b,
                      std::vector<double> &x, size_t max_iter, double tol);

int gmres(CSRMatrix &A, const std::vector<double> &b,
                      std::vector<double> &x, size_t max_iter, double tol);