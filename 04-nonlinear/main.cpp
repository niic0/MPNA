#include <iostream>
#include <stdint.h>
#include <mpi.h>

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

constexpr float k_0 = 0.01;
constexpr float sigma = 0.1; // 2. sigma = 1
constexpr float beta  = 1;   // 2. beta = 300
constexpr float gamma = 0.1; // Could be 1 or 10

int main(int argc, char* argv[]) {
  // MPI variable initialization
  MPI_INIT(&argc, &argv);

  int rank, num_procs;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_procs);

  // Partitionning the matrix by P mpi process
  int N = 100; // Matrix size
  int ilower = rank * (N / num_procs);
  int iupper = (rank + 1) * (N / num_procs) - 1;

  int jlower = ilower; // Square matrix
  int jupper = iupper; // Square matrix

  // Hypre matrix initilization
  HYPRE_IJMatrix ij_matrix;
  HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &ij_matrix);
  HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(ij_matrix);

  for(size_t i = 0; i<N; i++) 
  /* set matrix coefficients */
  //  HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, values);
  /* add-to matrix cofficients, if desired */
  //HYPRE_IJMatrixAddToValues(ij_matrix, nrows, ncols, rows, cols, values);

  //HYPRE_IJMatrixAssemble(ij_matrix);
  //HYPRE_IJMatrixGetObject(ij_matrix, (void **) &parcsr_matrix);


  // Gradiant conjugué non préconditionné

  /* solve */
  //  HYPRE_PCGSolve(solver, parcsr_A, b, x);

  HYPRE_Finalize(); /* must be the last HYPRE function call */

  return 0;
}
