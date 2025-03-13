#pragma once

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>

#include <iostream>
#include <mpi.h>
#include <vector>

/**
 * @brief Interface for handling HYPRE matrix assembly and RHS vector.
 */
class HypreInterface {
public:
  HypreInterface(size_t N, MPI_Comm comm);
  ~HypreInterface();

  void init();
  void assembleMatrix(const std::vector<std::vector<double>> &A);
  void assembleRHS(const std::vector<double> &F);
  void solveSystem(std::vector<double> &u);
  void finalize();

private:
  size_t N;       ///< Number of grid points
  MPI_Comm comm;  ///< MPI Communicator
  int rank, size; ///< MPI rank and total processes

  HYPRE_IJMatrix A_hypre;
  HYPRE_ParCSRMatrix parcsr_A;

  HYPRE_IJVector F_hypre;
  HYPRE_ParVector par_F;

  HYPRE_IJVector u_hypre;
  HYPRE_ParVector par_u;
};