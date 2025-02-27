#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Compressed Sparse Row (CSR) representation of a sparse matrix.
 */
class CSRMatrix {
private:
  unsigned int rows;
  unsigned int cols; ///< Number of columns
  int nnz;           ///< Number of non-zero elements

  size_t toMatId(size_t i, size_t j, size_t mesh_size) {
    return i + mesh_size * j;
  }

public:
  std::vector<double> values;        ///< Non-zero values (size = nnz)
  std::vector<size_t> columnIndices; ///< Column indices of non-zero values
  std::vector<size_t> rowOffsets;    ///< Row offsets (size = rows + 1)

  unsigned int getRows() { return rows; }
  unsigned int getCols() { return cols; }
  unsigned int getNnz() { return nnz; }

  /**
   * @brief Constructs a sparse matrix from a Matrix Market (.mtx) file.
   * @param filename The path to the .mtx file.
   */
  void fillWithFile(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Error: Unable to open file " + filename);
    }

    std::string line;

    // Skip comments (lines starting with '%')
    while (std::getline(file, line)) {
      if (line[0] != '%')
        break;
    }

    // Read matrix dimensions (rows, cols, nnz)
    std::istringstream iss(line);
    iss >> rows >> cols >> nnz;

    // Allocate fixed-size arrays
    values.reserve(nnz);
    columnIndices.reserve(nnz);
    rowOffsets.reserve(rows + 1);
    rowOffsets.push_back(0);

    // Temporary row counters to avoid if in the main loop
    std::vector<size_t> tempRowCount(rows, 0);

    size_t current_row = 0;
    size_t previous_row = -1; // -1 so the first iteration
                              // enter the if statement

    // Read matrix elements (row, col, value)
    for (size_t i = 0; i < nnz; i++) {
      int row, col;
      double value;
      file >> row >> col >> value;
      row--;
      col--; // Convert 1-based indexing to 0-based

      values.push_back(value);
      columnIndices.push_back(col);
      tempRowCount[row]++;
    }

    file.close();

    for (size_t i = 1; i < rows; i++) {
      rowOffsets.push_back(rowOffsets.back() + tempRowCount[i]);
    }
  }

  /**
   * @brief Fills a sparse matrix in CSR format for the Poisson equation using a
   * 5-point stencil.
   *
   * This function constructs the matrix representation of the discretized 2D
   * Poisson equation on a structured grid using finite differences. The matrix
   * is stored in Compressed Sparse Row (CSR) format.
   *
   * @param matrix The sparse matrix to be filled (in CSR format).
   * @param mesh_size The number of grid points along one dimension (total size
   * = mesh_size * mesh_size).
   */
  void fillWithPoisson(size_t mesh_size) {
    size_t n = mesh_size * mesh_size; // Total number of unknowns (matrix size)

    // Initialize matrix dimensions
    rows = n;
    cols = n;

    // Reserve memory to avoid multiple reallocations (5 non-zero
    // elements per row)
    values.reserve(5 * n - 4 * mesh_size);
    columnIndices.reserve(5 * n - 4 * mesh_size);
    rowOffsets.reserve(n + 1);

    size_t nnz = 0;          // Counter for non-zero elements
    rowOffsets.push_back(0); // First row starts at index 0

    // Loop over the 2D grid points
    for (size_t j = 0; j < mesh_size; j++) {
      for (size_t i = 0; i < mesh_size; i++) {
        size_t k = toMatId(i, j, mesh_size); // Convert (i, j) to matrix index

        // Diagonal element (central point)
        values.push_back(4.0);
        columnIndices.push_back(k);
        nnz++;

        // Left neighbor
        if (i > 0) {
          values.push_back(-1.0);
          columnIndices.push_back(toMatId(i - 1, j, mesh_size));
          nnz++;
        }

        // Right neighbor
        if (i < mesh_size - 1) {
          values.push_back(-1.0);
          columnIndices.push_back(toMatId(i + 1, j, mesh_size));
          nnz++;
        }

        // Bottom neighbor
        if (j > 0) {
          values.push_back(-1.0);
          columnIndices.push_back(toMatId(i, j - 1, mesh_size));
          nnz++;
        }

        // Top neighbor
        if (j < mesh_size - 1) {
          values.push_back(-1.0);
          columnIndices.push_back(toMatId(i, j + 1, mesh_size));
          nnz++;
        }

        // Mark the beginning of the next row
        rowOffsets.push_back(nnz);
      }
    }

    // Adjust the actual nnz value in case of boundary conditions
    nnz = nnz;
  }

  /**
   * @brief Computes the matrix-vector multiplication Ax = A * x using CSR
   * format.
   *
   * This function performs the sparse matrix-vector multiplication using the
   * Compressed Sparse Row (CSR) format.
   *
   * @param Ax The output vector (result of A * x).
   * @param x The input vector.
   */
  void matrixVectorMultiply(std::vector<double> &Ax,
                            const std::vector<double> &x) {
    size_t n = rows;

    // Ensure output vector Ax is correctly sized
    Ax.resize(n, 0.0);

    // Sparse matrix-vector multiplication
    for (size_t i = 0; i < n; i++) {
      double sum = 0.0;
      for (size_t j = rowOffsets[i]; j < rowOffsets[i + 1]; j++) {
        sum += values[j] * x[columnIndices[j]];
      }
      Ax[i] = sum;
    }
  }

  /**
   * @brief Prints the first rows and columns of a sparse matrix in CSR format.
   *
   * This function displays the matrix structure by iterating through its CSR
   * storage. It prints at most `max_rows` rows and `max_cols` columns.
   *
   * @param max_rows The maximum number of rows to print.
   * @param max_cols The maximum number of columns to print.x
   */
  void print(size_t max_rows, size_t max_cols) {
    size_t rows = std::min<size_t>(max_rows, rows); // Limit row display
    size_t cols = std::min<size_t>(max_cols, cols); // Limit column display

    std::cout << "Sparse Matrix (CSR) - Displaying first " << rows
              << " rows and " << cols << " columns:\n";

    for (size_t i = 0; i < rows; i++) {
      std::cout << "Row " << i << ": ";
      size_t row_start = rowOffsets[i];
      size_t row_end = rowOffsets[i + 1];

      size_t col_index = 0;
      for (size_t j = row_start; j < row_end; j++) {
        size_t col = columnIndices[j];

        // Only print values within the max_cols limit
        if (col < cols) {
          std::cout << "(" << col << ", " << values[j] << ") ";
        }
      }
      std::cout << "\n";
    }
  }

  /**
   * @brief Estimates the smallest eigenvalue of a sparse matrix using power
   * iteration.
   *
   * If the smallest eigenvalue is negative, the matrix is not positive
   * definite.
   *
   * @param A The sparse matrix in CSR format.
   * @return True if A is positive definite, false otherwise.
   */
  bool isPositiveDefinite() {
    size_t n = getRows();
    std::vector<double> x(n, 1.0), y(n, 0.0);

    for (int iter = 0; iter < 50; iter++) {
      matrixVectorMultiply(y, x);
      double lambda_min = *std::min_element(y.begin(), y.end());

      if (lambda_min < 0) {
        std::cout << "Matrix is not positive definite (smallest eigenvalue: "
                  << lambda_min << ")\n";
        return false;
      }

      // Normalize
      for (size_t i = 0; i < n; i++)
        x[i] = y[i] / lambda_min;
    }

    return true;
  }

  /**
   * @brief Checks if a sparse matrix in CSR format is symmetric.
   *
   * This function verifies that for every non-zero element A(i, j),
   * there exists a corresponding A(j, i) with the same value.
   *
   * @param A The sparse matrix in CSR format.
   * @return True if A is symmetric, false otherwise.
   */
  bool isSymmetric() {
    for (size_t i = 0; i < getRows(); i++) {
      for (size_t j = rowOffsets[i]; j < rowOffsets[i + 1]; j++) {
        size_t col = columnIndices[j];
        double value = values[j];
        bool found = false;

        for (size_t k = rowOffsets[col]; k < rowOffsets[col + 1]; k++) {
          if (columnIndices[k] == i && std::abs(values[k] - value) < 1e-9) {
            found = true;
            break;
          }
        }

        if (!found) {
          std::cout << "Matrix is not symmetric at (" << i << ", " << col
                    << ")\n";
          return false;
        }
      }
    }

    return true;
  }

  bool isDiagonalDominant() {
    for (size_t i = 0; i < getRows(); i++) {
      double diag = 0.0;
      double sum = 0.0;
      for (size_t j = rowOffsets[i]; j < rowOffsets[i + 1]; j++) {
        if (columnIndices[j] == i)
          diag = std::abs(values[j]);
        else
          sum += std::abs(values[j]);
      }
      if (diag < sum) {
        return false;
      }
    }

    return true;
  }
};
