#include "../include/CSRMatrix.h"
#include <cstdlib>
#include <mpi.h>
#include <vector>
#include <cmath>

/**
 * @brief Performs parallel matrix-vector multiplication using MPI with CSR
 * format.
 *
 * @param Ax The output vector (only valid on rank 0 after computation).
 * @param A The input matrix in CSR format (only fully stored on rank 0).
 * @param x The input vector.
 */
void matrixVectorMultiplyMPI(std::vector<double> &Ax, CSRMatrix &A,
                             const std::vector<double> &x) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t n;
  if (rank == 0) {
    n = A.getRows();
  }

  // Broadcast matrix size to all processes
  MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  // Calculate local row distribution
  size_t local_rows = n / size;
  size_t remainder = n % size;
  if (rank < remainder)
    local_rows++;

  size_t start_row =
      (rank < remainder) ? rank * local_rows : rank * local_rows + remainder;
  size_t end_row = start_row + local_rows;

  // Allocate storage on rank > 0
  if (rank != 0) {
    A.rowOffsets.resize(n + 1);
    A.columnIndices.resize(A.getNnz());
    A.values.resize(A.getNnz());
  }

  // Broadcast CSR matrix data
  MPI_Bcast(A.rowOffsets.data(), n + 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(A.columnIndices.data(), A.getNnz(), MPI_UNSIGNED_LONG, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(A.values.data(), A.getNnz(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Local computation
  std::vector<double> local_Ax(local_rows, 0.0);
  for (size_t i = start_row; i < end_row; i++) {
    double sum = 0.0;
    for (size_t j = A.rowOffsets[i]; j < A.rowOffsets[i + 1]; j++) {
      sum += A.values[j] * x[A.columnIndices[j]];
    }
    local_Ax[i - start_row] = sum;
  }

  // Gather results
  if (rank == 0) {
    Ax.resize(n, 0.0);
  }

  std::vector<int> recvcounts(size, local_rows);
  for (int i = 0; i < remainder; i++)
    recvcounts[i]++;

  std::vector<int> displs(size, 0);
  for (int i = 1; i < size; i++) {
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  }

  MPI_Gatherv(local_Ax.data(), local_rows, MPI_DOUBLE, Ax.data(),
              recvcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/**
 * @brief Computes the matrix-vector multiplication Ax = A * x using CSR
 * format.
 *
 * This function performs the sparse matrix-vector multiplication using the
 * Compressed Sparse Row (CSR) format.
 *
 * @param Ax The output vector (result of A * x).
 * @param A The input matrix.
 * @param x The input vector.
 */
void matrixVectorMultiply(std::vector<double> &Ax, CSRMatrix &A,
                          const std::vector<double> &x) {
  size_t n = A.getRows();

  // Ensure output vector Ax is correctly sized
  Ax.resize(n, 0.0);

  // Sparse matrix-vector multiplication
  for (size_t i = 0; i < n; i++) {
    double sum = 0.0;
    for (size_t j = A.rowOffsets[i]; j < A.rowOffsets[i + 1]; j++) {
      sum += A.values[j] * x[A.columnIndices[j]];
    }
    Ax[i] = sum;
  }
}

/**
 * @brief Creates a random filled vector
 *
 * @param v The Vector to fill
 * @param N The number of elements in the vector
 * @param seed The seed of the random generator
 */
void createVector(std::vector<double> &v, size_t N, int seed) {
  srand(seed);

  for (size_t i = 0; i < N; i++) {
    v[i] = (double)rand() / RAND_MAX;
  }
}

/**
 * @brief Computes the largest eigenvalue using the power iteration method (MPI
 * version).
 *
 * @param A The sparse matrix in CSR format.
 * @param max_iter Maximum number of iterations.
 * @param tol Convergence tolerance.
 * @return The estimated dominant eigenvalue.
 */
double powerIterationMPI(CSRMatrix &A, size_t max_iter, double tol) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t n = A.getRows();
  std::vector<double> x(n, 1.0); // Initial guess (all ones)
  std::vector<double> Ax(n, 0.0);

  double lambda_old = 0.0, lambda_new = 0.0;

  for (size_t k = 0; k < max_iter; k++) {
    // Perform matrix-vector multiplication Ax = A * x
    matrixVectorMultiplyMPI(Ax, A, x);

    // Compute the new estimated eigenvalue
    double local_dot = 0.0, local_norm = 0.0;
    for (size_t i = 0; i < n; i++) {
      local_dot += x[i] * Ax[i];   // x^T * Ax
      local_norm += Ax[i] * Ax[i]; // ||Ax||^2
    }

    // Sum across all MPI processes
    double global_dot, global_norm;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    lambda_new = global_dot / global_norm; // Rayleigh quotient

    // Normalize x
    double norm_factor = sqrt(global_norm);
    for (size_t i = 0; i < n; i++) {
      x[i] = Ax[i] / norm_factor;
    }

    // Check for convergence
    /*
    if (rank == 0) {
      std::cout << "Iteration " << k << ": λ = " << lambda_new << "\n";
    }*/

    if (std::abs(lambda_new - lambda_old) < tol) {
      break;
    }
    lambda_old = lambda_new;
  }

  return lambda_new;
}

/**
 * @brief Computes the Gershgorin bound for the eigenvalues of a matrix.
 *
 * This function estimates the largest possible eigenvalue using
 * Gershgorin's theorem. It computes the sum of absolute values of
 * non-diagonal elements for each row and determines the maximum bound.
 *
 * @param A The sparse matrix in CSR format.
 * @return The maximum Gershgorin bound.
 */
double computeGershgorinBound(CSRMatrix &A) {
  double max_bound = 0.0;
  for (size_t i = 0; i < A.getRows(); i++) {
    double center = 0.0, radius = 0.0;
    for (size_t j = A.rowOffsets[i]; j < A.rowOffsets[i + 1]; j++) {
      if (A.columnIndices[j] == i) {
        center = A.values[j]; // Diagonal element
      } else {
        radius += std::abs(A.values[j]); // Off-diagonal sum
      }
    }
    max_bound = std::max(max_bound, std::abs(center) + radius);
  }
  return max_bound;
}
