#include <iostream>
#include <iomanip>
#include <string>

/**
 * @brief Prints information about a sparse matrix before matrix-vector
 * multiplication.
 *
 * This function provides details about the matrix structure, such as its size,
 * number of non-zero elements, and sparsity ratio.
 *
 * @param filename Name of the matrix file (if applicable).
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @param nnz Number of non-zero elements.
 */
void printMatrixInfo(const std::string &filename, size_t rows, size_t cols, size_t nnz) {
    std::cout << std::endl << " 🔹 Matrix file               : " << filename << std::endl;
    std::cout << " 🔹 Matrix dimensions         : " << rows << " x " << cols << std::endl;
    std::cout << " 🔹 Number of non-zero values : " << nnz << std::endl;
    std::cout << " 🔹 Sparsity                  : " << std::fixed << std::setprecision(2)
              << (100.0 * nnz) / (rows * cols) << "%\n"
              << std::endl;
}

/**
 * @brief Prints performance results after a matrix-vector multiplication.
 *
 * @param time_ms Computation time in milliseconds.
 * @param mpi_nodes Number of MPI nodes used.
 */
void printMatrixVectorMultiplyMpiResults(double time_ms, size_t mpi_nodes) {
  std::cout << std::endl << "\e[1mMpi matrix-vector multiplication:\e[0m" << std::endl;

  std::cout << " MPI node   : " << mpi_nodes << std::endl;
  std::cout << " Time (ms)  : " << time_ms << std::endl;
}

/**
 * @brief Prints performance results after a matrix-vector multiplication.
 *
 * This function logs the execution time and MPI node count (if applicable).
 *
 * @param filename Name of the matrix file (if applicable).
 * @param time_ms Computation time in milliseconds.
 * @param mpi_nodes Number of MPI nodes used (default is 1 for non-MPI
 * execution).
 */
void printMatrixVectorMultiplySerialResults(double time_ms) {
  std::cout << std::endl << "\e[1mSerial matrix-vector multiplication:\e[0m" << std::endl;

  std::cout << " Time : " << time_ms << " ms" << std::endl;
}

/**
 * @brief Prints the results of the power iteration method.
 *
 * This function displays the dominant eigenvalue, computation time,
 * and the Gershgorin bound for validation.
 *
 * @param eigenvalue The computed dominant eigenvalue.
 * @param time_ms The execution time in milliseconds.
 * @param gershgorinBound The computed Gershgorin bound for eigenvalues.
 */
void printMatrixPowerIterationMpiResults(double eigenvalue, double time_ms,
                                         double gershgorinBound) {
  std::cout << std::endl << "\e[1mPower Iteration:\e[0m" << std::endl;

  std::cout << " Dominant Eigenvalue : " << eigenvalue << std::endl;
  std::cout << " Computation Time    : " << time_ms << " ms" << std::endl;
  std::cout << " Gershgorin Bound    : " << gershgorinBound << std::endl;
}