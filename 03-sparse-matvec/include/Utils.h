#pragma once

#include <cstddef>
#include <string>

void printMatrixInfo(const std::string &filename, size_t rows, size_t cols, size_t nnz);

void printMatrixVectorMultiplySerialResults(double time_ms);
void printMatrixVectorMultiplyMpiResults(double time_ms, size_t mpi_nodes);
void printMatrixPowerIterationMpiResults(double eigenvalue, double time_ms, double gershgorinBound);