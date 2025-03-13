#pragma once

#include <vector>
#include <cmath>
#include <iostream>

/**
 * @brief Builds the discretized nonlinear diffusion matrix and source term.
 */
class MatrixBuilder {
public:
    MatrixBuilder(size_t N, double dx, double k0, double sigma, double beta, double delta);

    void buildMatrix(std::vector<std::vector<double>> &A,
                                    const std::vector<double> &u);
    void buildRHS(std::vector<double> &F,
                                 const std::vector<double> &u);
    void printMatrix(const std::vector<std::vector<double>> &A) const;
    void printVector(const std::vector<double>& v) const;
    double kappa(double u) const;
    double Q(double x) const;

  private:
    size_t N;        ///< Number of grid points
    double dx;       ///< Grid spacing
    double k0;       ///< Conductivity coefficient
    double sigma;    ///< Radiation coefficient
    double beta;     ///< Source term coefficient
    double delta;    ///< Flame width
};
