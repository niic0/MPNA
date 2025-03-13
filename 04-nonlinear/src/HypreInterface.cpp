#include "../include/HypreInterface.h"

/**
 * @brief Constructor: Initializes HYPRE and MPI.
 */
HypreInterface::HypreInterface(size_t N, MPI_Comm comm) : N(N), comm(comm) {
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // HYPRE Initialization
  if (rank == 0) {
    HYPRE_Initialize();
  }

  // Compute local partitioning of matrix
  int local_size = N / size;
  int remainder = N % size;
  int ilower = rank * local_size + std::min(rank, remainder);
  int iupper = ilower + local_size - 1 + (rank < remainder);

  // Create HYPRE Matrix
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A_hypre);
  HYPRE_IJMatrixSetObjectType(A_hypre, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A_hypre);

  // Create HYPRE Vector for RHS
  HYPRE_IJVectorCreate(comm, ilower, iupper, &F_hypre);
  HYPRE_IJVectorSetObjectType(F_hypre, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(F_hypre);

  // Create HYPRE Vector for solution
  HYPRE_IJVectorCreate(comm, ilower, iupper, &u_hypre);
  HYPRE_IJVectorSetObjectType(u_hypre, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(u_hypre);
}

/**
 * @brief Destructor: Cleans up HYPRE objects.
 */
HypreInterface::~HypreInterface() {
  HYPRE_IJMatrixDestroy(A_hypre);
  HYPRE_IJVectorDestroy(F_hypre);
  HYPRE_IJVectorDestroy(u_hypre);
  HYPRE_Finalize();
}

/**
 * @brief Assemble matrix into HYPRE.
 */
void HypreInterface::assembleMatrix(const std::vector<std::vector<double>> &A) {
  int ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange(A_hypre, &ilower, &iupper, &jlower, &jupper);

  for (int i = ilower; i <= iupper; i++) {
    std::vector<int> cols;
    std::vector<double> values;

    for (size_t j = 0; j < N; j++) {
      if (A[i][j] != 0.0) {
        cols.push_back(j);
        values.push_back(A[i][j]);
      }
    }

    int nnz = cols.size();
    HYPRE_IJMatrixSetValues(A_hypre, 1, &nnz, &i, cols.data(), values.data());
  }

  HYPRE_IJMatrixAssemble(A_hypre);
  HYPRE_IJMatrixGetObject(A_hypre, (void **)&parcsr_A);
}

/**
 * @brief Assemble RHS vector into HYPRE.
 */
void HypreInterface::assembleRHS(const std::vector<double> &F) {
  int ilower, iupper;
  HYPRE_IJVectorGetLocalRange(F_hypre, &ilower, &iupper);

  std::vector<int> indices;
  std::vector<double> values;

  for (int i = ilower; i <= iupper; i++) {
    indices.push_back(i);
    values.push_back(F[i]);
  }

  HYPRE_IJVectorSetValues(F_hypre, values.size(), indices.data(),
                          values.data());
  HYPRE_IJVectorAssemble(F_hypre);
  HYPRE_IJVectorGetObject(F_hypre, (void **)&par_F);
}

/**
 * @brief Solve the linear system Ax = F using BoomerAMG.
 */
void HypreInterface::solveSystem(std::vector<double> &u) {
    HYPRE_Solver solver;
    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetTol(solver, 1e-7); // Tolerance
    HYPRE_BoomerAMGSetMaxIter(solver, 100); // Max iterations

    // 🛠️ Initialiser et remplir le vecteur `u_hypre` avec des valeurs initiales
    int ilower, iupper;
    HYPRE_IJVectorGetLocalRange(u_hypre, &ilower, &iupper);

    std::vector<int> indices(iupper - ilower + 1);
    std::vector<double> values(iupper - ilower + 1, 1.0); // Initial guess

    for (int i = ilower; i <= iupper; i++)
        indices[i - ilower] = i;
    
    // ⬇️ Mettre les valeurs initiales dans `u_hypre`
    HYPRE_IJVectorSetValues(u_hypre, values.size(), indices.data(), values.data());
    HYPRE_IJVectorAssemble(u_hypre);
    HYPRE_IJVectorGetObject(u_hypre, (void **) &par_u);  // 🔥 Nécessaire avant d'utiliser `par_u`

    // 🛠️ Lancer le solveur
    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_F, par_u);
    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_F, par_u);

    // 🔥 Récupérer la solution depuis `u_hypre`
    HYPRE_IJVectorGetValues(u_hypre, values.size(), indices.data(), values.data());

    // 🔥 Copier la solution dans `u`
    for (size_t i = 0; i < values.size(); i++)
        u[indices[i]] = values[i];

    // Nettoyage du solveur
    HYPRE_BoomerAMGDestroy(solver);
}


/**
 * @brief Finalizes HYPRE matrix and vector.
 */
void HypreInterface::finalize() {
  HYPRE_IJMatrixAssemble(A_hypre);
  HYPRE_IJVectorAssemble(F_hypre);
}
