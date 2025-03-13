# Nonlinear Diffusion Solver Using HYPRE
## Overview
This project implements a nonlinear diffusion solver for temperature distribution in a flame using finite difference discretization and implicit numerical schemes. The solver utilizes the **HYPRE** library. The implementation supports both **linearized implicit schemes** and **Newton-Raphson iterations** for solving the nonlinear equations.

---

## Build Instructions
The project uses **CMake** and **Ninja** for compilation.

### Dependencies
- **HYPRE** ([https://github.com/hypre-space/hypre](https://github.com/hypre-space/hypre))
- **MPI**
- **CMake**
- **Ninja**

### Building the Project
Create a build directory:
```bash
mkdir build && cd build
```

Run CMake with Ninja:
```bash
cmake -G Ninja ..
```

Compile the project:
```bash
ninja
```

Run the executable with MPI:
```bash
mpirun -np <num_processes> ./diffusion
```

## Usage and Parameters
The solver can be customized through several parameters defined in the source code.

### User Parameters
These parameters can be modified in **main.cpp**:

| Parameter   | Description                        |
| ----------- | ---------------------------------- |
| `N`         | Number of grid points              |
| `L`         | Domain length (default = 1.0)      |
| `dx`        | Grid spacing (`dx = L/(N-1)`)      |
| `k0`        | Thermal conductivity constant      |
| `sigma`     | Radiation parameter                |
| `beta`      | Heat source parameter              |
| `gamma`     | Time-stepping control (0.1, 1, 10) |
| `tol`       | Convergence tolerance              |
| `max_iters` | Maximum number of iterations       |

### Output
The program prints the iteration details to the console:

```
Iteration 0  u_max: 1.1  dt: 0.0017132  Error: 0.0882982
Iteration 1  u_max: 1.72077  dt: 0.000698282  Error: 0.384564
...
```

Where:
- `u_max` is the maximum temperature at each step.
- `dt` is the computed time step.
- `Error` is the difference from the previous iteration.

## Mathematical Formulation
The governing equation for nonlinear diffusion is given by:

$$
-frac{\partial}{\partial x} \left( \kappa(u) frac{\partial u}{\partial x}) + \right (sigma (u^4 - 1)) = Q(x)
$$

where:

- $\kappa(u) = \kappa_0 u^q$ is the temperature-dependent thermal conductivity.
- $Q(x)$ represents the heat source (a step function).

### Finite Difference Discretization
The PDE is discretized using central differences:

$$
\left(frac{\kappa_{i+frac{1}{2}} (u_{i+1} - u_i) - \kappa_{i-frac{1}{2}} (u_i - u_{i-1})}{dx^2}) + \sigma (u_i^4 - 1) = Q_i
$$

Boundary conditions are applied:
- *Neumann* at $x = 0$: $\frac{du}{dx}(0) = 0$
- *Dirichlet* at $x = 1$: $u(1) = 1$

### Implicit Scheme Implementation
The implicit scheme solves:

$$
A u^{n+1} = u^n + dt (Q_i + \sigma)
$$

where $A$ is the system matrix stored in **HYPRE** using the IJ interface. The system is solved using:

- **BoomerAMG** (Algebraic Multigrid Preconditioner)
- **PCG** (Preconditioned Conjugate Gradient)

### Newton-Raphson Method
For full nonlinear solving, Newton-Raphson iterations are used:

1. Compute the Jacobian matrix:
   $J_{i,j} = \frac{\partial F_i}{\partial u_j}$
2. Solve the Newton update equation:
   $u^{k+1} = u^k - J^{-1} F(u^k)$
3. Use **HYPRE** to solve the Jacobian system.

Newton’s method is compared against the implicit scheme for convergence speed and residual decay.

---

## HYPRE Integration
HYPRE is used for efficiently solving large sparse linear systems. This project employs:

- **IJ Interface**: Handles matrix storage.
- **BoomerAMG**: Multigrid preconditioner for fast convergence.

### HYPRE Implementation

1. Assemble the matrix in HYPRE:
   ```cpp
   hypre.assembleMatrix(A);
   hypre.assembleRHS(F);
   hypre.finalize();
   ```
2. Solve the system using BoomerAMG:
   ```cpp
   hypre.solveSystem(u);
   ```
3. Retrieve the solution:
   ```cpp
   MPI_Allreduce(u.data(), global_u.data(), N, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   u = global_u;
   ```