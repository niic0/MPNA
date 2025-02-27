#include <cmath>

#include "../include/CSRMatrix.h"


/**
 * @brief Checks if a matrix is suitable for Conjugate Gradient (CG).
 *
 * @param A The sparse matrix in CSR format.
 * @return True if CG can be used, false otherwise.
 */
bool canUseConjugateGradient(CSRMatrix &A) {
    return A.isSymmetric() && A.isPositiveDefinite();
}

/**
 * @brief Solves the linear system Ax = b using the Conjugate Gradient (CG)
 * method.
 *
 * This function applies the CG algorithm to solve a symmetric positive-definite
 * linear system stored in CSR format. It iterates until the residual is below
 * a given tolerance or the maximum number of iterations is reached.
 *
 * @param A The sparse matrix stored in CSR format.
 * @param b The right-hand side vector.
 * @param x The solution vector (updated in-place).
 * @param max_iter Maximum number of iterations before stopping.
 * @param tol Convergence tolerance.
 *
 * @return Return the number of iteration if converged, -1 if not.
 */
int conjugateGradient(CSRMatrix &A, const std::vector<double> &b,
                      std::vector<double> &x, size_t max_iter, double tol) {
  if (canUseConjugateGradient(A) == false) {
    std::cout << "Can't use Conjugate Gradiant on this matrix !" << std::endl;
    return -1;
  }

  size_t n = A.getRows();

  // Initialize residual, direction, and temporary vectors
  std::vector<double> r(n, 0.0);
  std::vector<double> p(n, 0.0);
  std::vector<double> Ap(n, 0.0);

  // Compute initial residual r0 = b - A * x
  A.matrixVectorMultiply(Ap, x);
  for (size_t i = 0; i < n; i++) {
    r[i] = b[i] - Ap[i];
    p[i] = r[i]; // Initial direction p0 = r0
  }

  double r_norm_old = 0.0;
  for (size_t i = 0; i < n; i++)
    r_norm_old += r[i] * r[i]; // r^T * r

  for (unsigned int k = 0; k < max_iter; k++) {
    // Compute Ap = A * p
    A.matrixVectorMultiply(Ap, p);

    // Compute alpha_k = (r_k^T * r_k) / (p_k^T * A * p_k)
    double pAp = 0.0;
    for (size_t i = 0; i < n; i++)
      pAp += p[i] * Ap[i];

    double alpha = r_norm_old / pAp;

    // Update solution x_k+1 = x_k + alpha_k * p_k
    for (size_t i = 0; i < n; i++)
      x[i] += alpha * p[i];

    // Update residual r_k+1 = r_k - alpha_k * A * p_k
    double r_norm_new = 0.0;
    for (size_t i = 0; i < n; i++) {
      r[i] -= alpha * Ap[i];
      r_norm_new += r[i] * r[i]; // Compute ||r_k+1||^2
    }

    // Check convergence
    if (std::sqrt(r_norm_new) < tol) {
      return k+1;
    }

    // Compute beta_k = (r_k+1^T * r_k+1) / (r_k^T * r_k)
    double beta = r_norm_new / r_norm_old;
    r_norm_old = r_norm_new;

    // Update direction p_k+1 = r_k+1 + beta_k * p_k
    for (size_t i = 0; i < n; i++)
      p[i] = r[i] + beta * p[i];
  }

  return -1;
}