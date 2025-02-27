# Sparse Linear Solvers

## Project Overview
This project implements **iterative solvers** (Jacobi, Gauss-Seidel, GMRES, Conjugate Gradient) to solve **sparse linear systems** of the form:

$$
Ax = b
$$

where $A$ is a sparse matrix stored in **Compressed Sparse Row (CSR) format**.

---

## Methodology

### Sparse Matrix Representation (CSR)
To efficiently store sparse matrices, this project uses the **Compressed Sparse Row (CSR) format**, which consists of:
- **values**: Stores the non-zero elements.
- **columnIndices**: Stores the column indices of non-zero elements.
- **rowOffsets**: Stores the starting position of each row in the `values` array.

### Iterative Solvers Implemented
This project includes **four iterative solvers**:

1. **Jacobi Method**  

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)
$$

   - Converges for **diagonally dominant** matrices.
   - Requires **global communication** in MPI.

2. **Gauss-Seidel Method**  

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
$$

   - Faster than Jacobi, as updates are **immediately used**.
   - Converges for **symmetric positive definite (SPD) matrices**.

3. **Conjugate Gradient (CG) Method**  
   - Best suited for **SPD matrices**.
   - Uses **search directions** instead of simple updates:

     $$
     p_k = r_k + \beta_k p_{k-1}
     $$

   - Requires **matrix-vector products** in CSR format.

4. **Generalized Minimal Residual (GMRES) Method**  
   - Works for **non-symmetric matrices**.
   - Constructs an **orthonormal basis** using **Arnoldi iteration**.
   - More expensive than CG but handles a **wider range of problems**.

---

## Installation & Compilation
This project requires:
- **CMake** for compilation

### Compile the Project
```bash
mkdir build && cd build
cmake .. -G Ninja
cmake --build .
```

## Usage
```bash
Usage: ./SparseSolver [OPTIONS]
Options:
  -i, --input <file>       Path to input matrix file (.mtx)
  -m, --mesh <size>        Generate a Poisson matrix with given mesh size
  -n, --iterations <num>   Set maximum iterations (default: 100000)
  -t, --tolerance <value>  Set convergence tolerance (default: 1e-6)
  -h, --help               Display this help message
```

Examples:
```bash
# Run solver on Poisson generated matrix of mesh size 50
./SparseSolvers -m 50
```

## Project Structure
```bash
├── build/
├── include/
│   ├── ArgumentParser.h     # Parses user arguments
│   ├── Benchmark.h
│   ├── CSRMatrix.h          # Sparse matrix class
│   ├── Solver.h
│   ├── Utils.h
├── input/                   # Matrix market input files
│   ├── bcsstk03.mtx
│   ├── cfd1.mtx
│   ├── example.mtx
├── src/
│   ├── ArgumentParser.cpp
│   ├── Benchmark.cpp
│   ├── ConjugateGradient.cpp
│   ├── GaussSeidel.cpp
│   ├── GMRES.cpp
│   ├── Jacobi.cpp
│   ├── Utils.cpp
├── main.cpp                 # Main entry point
├── README.md
```
