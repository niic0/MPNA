# MPI-Based Sparse Matrix Power Iteration Solver

## 📌 Project Overview
This project implements **Power Iteration** to compute the **dominant eigenvalue** of a sparse matrix using **MPI**. The matrix is stored in **Compressed Sparse Row (CSR) format**.

---

## Installation & Compilation
This project requires:
- **MPI** for distributed computation
- **CMake** for compilation

### Compile the Project
```bash
mkdir build && cd build
cmake .. -G Ninja
cmake --build .
```

## Usage
1. Running the Serial Version
```bash
./matvec -i ../input/bcsstk03.mtx
```

2. Running the MPI Version
Run with different number of processes:
```bash
mpirun -np 1 ./matvec -i ../input/bcsstk03.mtx
mpirun -np 2 ./matvec -i ../input/bcsstk03.mtx
mpirun -np 4 ./matvec -i ../input/bcsstk03.mtx
```

## Methodology

### Power Iteration
Power Iteration finds the **largest magnitude eigenvalue** $\lambda_{\max}$ of a matrix $A$ by iterating:

$$
x_{k+1} = \frac{Ax_k}{\|Ax_k\|}
$$

The estimated eigenvalue is computed using the **Rayleigh quotient**:
$$
\lambda_k = \frac{x_k^T A x_k}{x_k^T x_k}
$$

#### Convergence Criteria
The algorithm stops when:
$$
|\lambda_{k+1} - \lambda_k| < \text{tolerance}
$$
where the tolerance is user-defined.

### Correctness Check: Gershgorin’s Theorem
To verify the computed eigenvalue, we compare it with the **Gershgorin bound**, which states that **all eigenvalues** of $A$ lie within disks centered at $a_{ii}$ with radius:
$$
R_i = \sum_{j \neq i} |a_{ij}|
$$

If $\lambda_{\max}$ is within this bound, the result is valid.

---

## CSR Matrix Format
To efficiently store sparse matrices, this project uses the **Compressed Sparse Row (CSR) format**, which consists of:
- **values**: Non-zero elements of the matrix.
- **columnIndices**: Column indices of the non-zero elements.
- **rowOffsets**: Starting index of each row in the `values` array.

## Project Structure
```
├── build/
├── include/
│   ├── CSRMatrix.h          # Sparse matrix class
│   ├── MatrixVector.h
│   ├── Utils.h
├── input/                   # Matrix market input files
├── src/
│   ├── ConjugateGradient.cpp
│   ├── GaussSeidel.cpp
│   ├── Jacobi.cpp
│   ├── MatrixVector.cpp
│   ├── Utils.cpp
├── main.cpp                 # Main entry point
├── README.md
```