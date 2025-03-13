#include "../include/MatrixBuilder.h"

#include "../include/MatrixBuilder.h"

MatrixBuilder::MatrixBuilder(size_t N, double dx, double k0, double sigma,
                             double beta, double delta)
    : N(N), dx(dx), k0(k0), sigma(sigma), beta(beta), delta(delta) {}

/**
 * @brief Computes the diffusion coefficient \(\kappa(u)\).
 */
double MatrixBuilder::kappa(double u) const {
  return k0 * std::pow(u, 2); // Example choice: κ(u) = k0 * u²
}

/**
 * @brief Computes the source term \( Q(x) \) using a Heaviside function.
 */
double MatrixBuilder::Q(double x) const { return (x < delta) ? beta : 0.0; }

/**
 * @brief Builds the coefficient matrix A for the nonlinear diffusion problem.
 *
 * This function fills the tridiagonal matrix A using finite difference
 * discretization.
 *
 * @param A Reference to the matrix to be constructed.
 * @param u Reference to the temperature field.
 */
void MatrixBuilder::buildMatrix(std::vector<std::vector<double>> &A,
                                const std::vector<double> &u) {
  A.assign(N, std::vector<double>(N, 0.0)); // Initialize NxN matrix with zeros
  double dx2 = dx * dx;

  for (size_t i = 1; i < N - 1; i++) {
    // Compute diffusion coefficients at i±1/2
    double kappa_ip = (kappa(u[i + 1]) + kappa(u[i])) / 2.0;
    double kappa_im = (kappa(u[i - 1]) + kappa(u[i])) / 2.0;

    // Construct the matrix A
    A[i][i - 1] = -kappa_im / dx2;
    A[i][i] = (kappa_im + kappa_ip) / dx2 + 4.0 * sigma * std::pow(u[i], 3);
    A[i][i + 1] = -kappa_ip / dx2;
  }

  // Apply boundary conditions
  A[0][0] = 1.0; // Neumann (mirroring)
  A[0][1] = -1.0;

  A[N - 1][N - 1] = 1.0; // Dirichlet
}

/**
 * @brief Builds the right-hand side vector F(u).
 *
 * This function fills the RHS vector F based on the nonlinear source term
 * and radiation effect.
 *
 * @param F Reference to the right-hand side vector.
 * @param u Reference to the temperature field.
 */
void MatrixBuilder::buildRHS(std::vector<double> &F,
                             const std::vector<double> &u) {
  F.assign(N, 0.0); // Initialize RHS vector with zeros

  for (size_t i = 1; i < N - 1; i++) {
    double x = i * dx;
    F[i] = sigma * (std::pow(u[i], 4) - 1.0) + Q(x);
  }

  // Apply boundary conditions
  F[0] = 0.0;     // Neumann (mirroring)
  F[N - 1] = 1.0; // Dirichlet
}

/**
 * @brief Prints the matrix.
 */
void MatrixBuilder::printMatrix(
    const std::vector<std::vector<double>> &A) const {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }
}

/**
 * @brief Prints a vector.
 */
void MatrixBuilder::printVector(const std::vector<double> &v) const {
  for (size_t i = 0; i < v.size(); i++) {
    std::cout << v[i] << " ";
  }
  std::cout << "\n";
}
