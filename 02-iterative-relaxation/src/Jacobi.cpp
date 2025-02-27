#include <vector>

#include "../include/CSRMatrix.h"

/**
 * @brief
 *
 */
int jacobi(CSRMatrix &A, const std::vector<double> &b, std::vector<double> &x,
           size_t max_iter, double tol) {
  // Verify the matrix diaognal dominance
  if(A.isDiagonalDominant() == false){
    std::cout << "Can't use Jacobi on non diognal dominant matrix !"
              << std::endl;
  }

  std::vector<double> x_new(A.getRows(), 0.0);

  for (size_t k = 0; k < max_iter; k++) {
    double maxDiff = 0.0;

    // Sum aii * xj^(k)
    for (size_t i = 0; i < A.getRows(); i++) {
      double sum = 0.0;
      double diag_val = 0.0;

      for (size_t j = A.rowOffsets[i]; j < A.rowOffsets[i + 1]; j++) {
        int col = A.columnIndices[j];
        double val = A.values[j];

        if (col == i) {
          diag_val = val;
        } else {
          sum += val * x[col];
        }
      }

      // update xi^(k+1)
      x_new[i] = (b[i] - sum) / diag_val;

      // Verify convergence
      maxDiff = std::max(maxDiff, std::abs(x_new[i] - x[i]));
    }

    // update x vector
    x = x_new;

    // Verify convergence with the tolerance and return the number of iteration
    if (maxDiff < tol) {
      return k + 1;
    }
  }
  
  return -1;
}