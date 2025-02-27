# MPNA Courses

## Overview
This repository contains **practical assignments** for the course **Méthodes de Programmation Numérique Avancées (MPNA)**. Each assignment focuses on a different aspect of **numerical algorithms** for high-performance computing, including matrix operations, iterative solvers, sparse matrix computations, and nonlinear systems.

## Repository Structure
The repository is organized into four main practical assignments:
```bash
├── 01-matrix-matrix        # Matrix-Matrix Operations 
├── 02-iterative-relaxation # Iterative Methods (Jacobi, Gauss-Seidel, GMRES, CG) 
├── 03-sparse-matvec        # Sparse Matrix-Vector Multiplication 
└── 04-nonlinear            # Nonlinear Equation Solvers using HYPRE
```

Each folder contains:
- **Source Code** (only `c++`)
- **Makefiles / CMakeLists.txt / Ninja** for compilation
- **Input Data** (Matrix Market files)
- **Documentation** (`README.md` per assignment)