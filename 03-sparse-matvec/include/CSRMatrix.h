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
   * @brief Resizes the CSR matrix to allocate memory for MPI.
   *
   * This method is used on non-root ranks to allocate space for incoming data.
   *
   * @param new_rows Number of rows.
   * @param new_cols Number of columns.
   * @param new_nnz  Number of non-zero values.
   */
  void resize(size_t new_rows, size_t new_cols, size_t new_nnz) {
    rows = new_rows;
    cols = new_cols;
    nnz = new_nnz;

    // Allocate space for CSR structure
    values.resize(nnz);
    columnIndices.resize(nnz);
    rowOffsets.resize(rows + 1);
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
};
