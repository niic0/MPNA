

#include "../include/CSRMatrix.h"

/**
 * @brief Solves a linear system Ax = b using the Gauss-Seidel iterative method
 * with a CSR matrix.
 *
 * This function applies the Gauss-Seidel algorithm, which updates the solution
 * vector x in-place by iterating over each row and solving for x_i using the
 * latest available values.
 *
 * @param A The sparse matrix stored in Compressed Sparse Row (CSR) format.
 * @param b The right-hand side vector.
 * @param x The solution vector (updated in-place).
 * @param max_iters Maximum number of iterations before stopping.
 * @param tol Convergence tolerance (stops when the solution change is below
 * this value).
 *
 * @return Return the number of iteration if converged, -1 if not.
 */
int gaussSeidel(CSRMatrix &A, const std::vector<double> &b,
                std::vector<double> &x, size_t max_iter, double tol) {
  size_t n = A.getRows(); // Number of rows
  for (unsigned int iter = 0; iter < max_iter; ++iter) {
    bool converged = true;

    // Iterate over each row
    for (size_t i = 0; i < n; ++i) {
      double sum = b[i]; // Initialize with b[i]
      double diag = 0.0; // Stores the diagonal value

      // Iterate over non-zero elements in row i
      for (size_t idx = A.rowOffsets[i]; idx < A.rowOffsets[i + 1]; ++idx) {
        size_t col = A.columnIndices[idx];
        double value = A.values[idx];

        if (col == i) {
          diag = value; // Extract diagonal value
        } else {
          sum -= value * x[col]; // Compute the sum using previous values
        }
      }

      // Compute new value of x[i]
      double new_x_i = sum / diag;

      // Check for convergence
      if (std::abs(new_x_i - x[i]) > tol) {
        converged = false;
      }

      x[i] = new_x_i; // Update x[i]
    }

    // Stop if the solution has converged
    if (converged) {
      return iter + 1;
    }
  }

  return -1;
}