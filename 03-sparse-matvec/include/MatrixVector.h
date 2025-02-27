#pragma once

#include "CSRMatrix.h"

void matrixVectorMultiplyMPI(std::vector<double> &Ax, CSRMatrix &A,
                             const std::vector<double> &x);
void matrixVectorMultiply(std::vector<double> &Ax, CSRMatrix &A,
                          const std::vector<double> &x);

void createVector(std::vector<double> &v, size_t N, int seed);

double powerIterationMPI(CSRMatrix &A, size_t max_iter, double tol);

double computeGershgorinBound(CSRMatrix &A);