#include "../include/CSRMatrix.h"
#include <cmath>
#include <vector>

/**
 * @brief Computes the dot product of two vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The dot product.
 */
static double dot_product(const std::vector<double> &v1,
                          const std::vector<double> &v2) {
  double result = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

/**
 * @brief Computes the Euclidean norm of a vector.
 *
 * @param v The vector.
 * @return The Euclidean norm.
 */
static double norm(const std::vector<double> &v) {
  return std::sqrt(dot_product(v, v));
}

/**
 * @brief Solves a linear system Ax = b using the restarted GMRES method.
 *
 * This function implements the GMRES algorithm using an Arnoldi process with
 * Givens rotations to solve a non-symmetric linear system stored in CSR format.
 *
 * @param A The sparse matrix stored in CSR format.
 * @param b The right-hand side vector.
 * @param x The initial guess and eventual solution vector (updated in-place).
 * @param max_iter Maximum number of iterations allowed.
 * @param tol Convergence tolerance.
 * @return The total number of iterations if converged, -1 otherwise.
 */
int gmres(CSRMatrix &A, const std::vector<double> &b, std::vector<double> &x,
          size_t max_iter, double tol) {
  const size_t n = A.getRows();
  // Use a fixed restart value (e.g., 50) or max_iter if it is lower.
  const size_t restart = (max_iter < 50) ? max_iter : 50;
  size_t total_iters = 0;

  // Temporary vector to hold A*x
  std::vector<double> Ax(n, 0.0);

  // Compute initial residual r = b - A*x
  A.matrixVectorMultiply(Ax, x);
  std::vector<double> r(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    r[i] = b[i] - Ax[i];
  }
  double beta = norm(r);
  if (beta < tol) {
    return total_iters;
  }

  // Restart loop
  while (total_iters < max_iter) {
    // Recompute residual r = b - A*x and its norm
    A.matrixVectorMultiply(Ax, x);
    for (size_t i = 0; i < n; ++i) {
      r[i] = b[i] - Ax[i];
    }
    beta = norm(r);
    if (beta < tol) {
      return total_iters;
    }

    // Allocate Krylov basis V with (restart+1) vectors
    std::vector<std::vector<double>> V(restart + 1,
                                       std::vector<double>(n, 0.0));
    // Normalize r to obtain the first basis vector V[0]
    for (size_t i = 0; i < n; ++i) {
      V[0][i] = r[i] / beta;
    }

    // Allocate the Hessenberg matrix H with dimensions (restart+1) x restart
    std::vector<std::vector<double>> H(restart + 1,
                                       std::vector<double>(restart, 0.0));

    // Vectors for storing Givens rotation coefficients
    std::vector<double> cs(restart, 0.0);
    std::vector<double> sn(restart, 0.0);

    // g vector holds the residual norm projections
    std::vector<double> g(restart + 1, 0.0);
    g[0] = beta;

    size_t j = 0;
    // Arnoldi process with Givens rotations
    for (; j < restart && total_iters + j < max_iter; ++j) {
      // Compute w = A * V[j]
      std::vector<double> w(n, 0.0);
      A.matrixVectorMultiply(w, V[j]);

      // Orthogonalize w against previous basis vectors V[0] ... V[j]
      for (size_t i = 0; i <= j; ++i) {
        H[i][j] = dot_product(w, V[i]);
        // Subtract the projection from w
        for (size_t k = 0; k < n; ++k) {
          w[k] -= H[i][j] * V[i][k];
        }
      }

      // Compute the norm of w (this becomes H[j+1][j])
      H[j + 1][j] = norm(w);
      // Normalize w to form the next basis vector
      if (H[j + 1][j] != 0.0) {
        for (size_t i = 0; i < n; ++i) {
          V[j + 1][i] = w[i] / H[j + 1][j];
        }
      }

      // Apply all previous Givens rotations to the new column H[0..j+1][j]
      for (size_t i = 0; i < j; ++i) {
        double temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
        H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
        H[i][j] = temp;
      }

      // Compute the new Givens rotation to eliminate H[j+1][j]
      double delta = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
      if (delta == 0.0) {
        cs[j] = 1.0;
        sn[j] = 0.0;
      } else {
        cs[j] = H[j][j] / delta;
        sn[j] = H[j + 1][j] / delta;
      }
      // Apply the rotation to H[j][j] and update the residual vector g
      H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
      g[j + 1] = -sn[j] * g[j];
      g[j] = cs[j] * g[j];

      // Check convergence using the updated residual
      if (std::abs(g[j + 1]) < tol) {
        ++j; // Count the current iteration
        break;
      }
    }

    // Solve the least-squares problem H*y = g via back substitution.
    // The system is upper triangular with j equations.
    std::vector<double> y(j, 0.0);
    for (int i = static_cast<int>(j) - 1; i >= 0; --i) {
      y[i] = g[i];
      for (size_t k = i + 1; k < j; ++k) {
        y[i] -= H[i][k] * y[k];
      }
      y[i] /= H[i][i];
    }

    // Update the solution: x = x + V[0:j]*y
    for (size_t i = 0; i < n; ++i) {
      double update = 0.0;
      for (size_t k = 0; k < j; ++k) {
        update += y[k] * V[k][i];
      }
      x[i] += update;
    }

    total_iters += j;

    // Check final residual norm for convergence
    A.matrixVectorMultiply(Ax, x);
    for (size_t i = 0; i < n; ++i) {
      r[i] = b[i] - Ax[i];
    }
    beta = norm(r);
    if (beta < tol) {
      return total_iters;
    }
  }

  // Return -1 if the solver did not converge within max_iter iterations.
  return -1;
}
