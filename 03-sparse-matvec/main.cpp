#include "include/ArgumentParser.h"
#include "include/CSRMatrix.h"
#include "include/MatrixVector.h"
#include "include/Utils.h"
#include <chrono>
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CliConfig config;
  parseArguments(argc, argv, config);

  // Load matrix on rank 0
  CSRMatrix A;
  if (rank == 0) {
    A.fillWithFile(config.input_file);
  }

  // Broadcast matrix structure
  size_t rows, cols, nnz;
  if (rank == 0) {
    rows = A.getRows();
    cols = A.getCols();
    nnz = A.getNnz();
  }
  MPI_Bcast(&rows, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cols, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nnz, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  // Resize on non-root ranks
  if (rank != 0) {
    A.resize(rows, cols, nnz);
  }

  // Broadcast full CSR data
  MPI_Bcast(A.rowOffsets.data(), rows + 1, MPI_UNSIGNED_LONG, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(A.columnIndices.data(), nnz, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(A.values.data(), nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<double> x(rows, 1.0);
  std::vector<double> Ax_serial, Ax_parallel;

  // Run serial matrix vector multiplication
  if (rank == 0) {
    printMatrixInfo(config.input_file, rows, cols, nnz);

    auto start = std::chrono::high_resolution_clock::now();
    matrixVectorMultiply(Ax_serial, A, x);
    auto end = std::chrono::high_resolution_clock::now();
    double time_serial =
        std::chrono::duration<double, std::milli>(end - start).count();

    printMatrixVectorMultiplySerialResults(time_serial);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Run Mpi matrix vector multiplication
  auto start = std::chrono::high_resolution_clock::now();
  matrixVectorMultiplyMPI(Ax_parallel, A, x);
  auto end = std::chrono::high_resolution_clock::now();
  double time_parallel =
      std::chrono::duration<double, std::milli>(end - start).count();

  if (rank == 0) {
    printMatrixVectorMultiplyMpiResults(time_parallel, size);
  }

  // Run Power Iteration
  start = std::chrono::high_resolution_clock::now();
  double eigenvalue = powerIterationMPI(A, 1000, 1e-6);
  end = std::chrono::high_resolution_clock::now();
  double time_ms =
      std::chrono::duration<double, std::milli>(end - start).count();

  if (rank == 0) {
    double gershgorinBound = computeGershgorinBound(A);
    printMatrixPowerIterationMpiResults(eigenvalue, time_ms, gershgorinBound);
  }

  MPI_Finalize();

  return 0;
}
