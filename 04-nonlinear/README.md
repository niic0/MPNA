# Nonlinear Diffusion Equation

Source: https://perso.univ-lyon1.fr/marc.buffat/COURS/COURSDF_HTML/node29.html

## Physical Problem: Diffusion in a Flame

We want to calculate the transverse temperature distribution in a flame. For this, we use a simple model, assuming that in the transverse direction (denoted $x$) there are heat exchanges only by diffusion and radiation.

![flamme](flamme.png)

With these assumptions, the stationary energy conservation equation in the $x$ direction is written as:

$$ \underbrace{-\frac{\partial}{\partial}\left(\lambda\frac{\partial T}{\partial x}\right)}_{\mbox{diffusion}} + \underbrace{\sigma(T^{4}-T_{0}^{4})}_{\mbox{rayonnement}}=\underbrace{Q}_{\mbox{source}} $$

  
In this model, the energy $Q$ produced in the flame of width $2\delta$ by chemical reaction is diffused by conduction and radiated to the outside air at temperature $T_{0}$. For radiation, we have adopted a simple black body radiation model proportional to $T^{4}$ with a radiation constant $\sigma$. The temperature variations between the flame and the outside are significant (with a ratio of about 3 to 4), the conduction coefficient $\lambda$ depends on the temperature through a power law $\lambda(T)=\lambda_{0}T^{q}$. Finally, for the source term, we will choose a constant reaction term in the flame and zero outside.
$$\frac{\partial T}{\partial x}(x=0)=0$$ 
and
$$T(x=L)=T_{0}.$$ 

$L$ étant une distance grande par rapport à l'épaisseur $\delta$ de la flamme.

En posant $u=\frac{T}{T_{0}}$, le problème modèle associé s'écrit alors (en choisissant $L=1$ et $\delta=0.2$):

$$-\frac{\partial}{\partial}\left(\kappa(u)\frac{\partial u}{\partial x}\right)+\sigma(u^{4}-1)=Q(x),       x\in]0,1[$$

$$\frac{\partial u}{\partial x}(0)=0, u(1)=1$$

with $$Q(x)=\beta  Heaviside(0.2-x)$$ and $$\kappa(u)=\kappa_{0}u^{q}.$$


### Discretization by Finite Differences

To numerically solve the nonlinear problem, we first discretize it using centered finite differences with a mesh size of $dx=1/N$ in the form:

$$-\left(\frac{\kappa_{i+\frac{1}{2}}\left(u_{i+1}-u_{i}\right)-\kappa_{i-\frac{1}{2}}\left(u_{i}-u_{i-1}\right)}{dx^{2}}\right)+\sigma\left(u_{i}^{4}-1\right)=Q_{i}    \forall i=0,N-1$$

noting that $\kappa_{i+\frac{1}{2}}=\kappa\frac{\kappa(u_{i+1})+\kappa(u_{i})}{2}$ and $\kappa_{i-\frac{1}{2}}=\frac{\kappa(u_{i-1})+\kappa(u_{i})}{2}$ are the diffusion coefficients at $i+\frac{1}{2}$ and $i-\frac{1}{2}$.

To these equations, we add the 2 boundary conditions:

1. the Neumann condition at $i=0$, translated as a mirror condition for the equation at $i=0$: $u_{-1}=u_{1}$,
2. the Dirichlet condition at $i=N$: $u_{N}=1$.

Thus, we obtain $N+1$ nonlinear equations for $N+1$ unknowns $\{u_{i}\}_{i=0,N}$, which we symbolically write as:
$$F_{i}(u_{0},u_{1},..u_{j}...,u_{N})=0, \forall i=0,N$$

To numerically solve these equations, we transform the problem into an equivalent fixed-point problem:

$$u_{i}=G_{i}(u_{0},u_{1},..u_{j}...,u_{N}), \forall i=0,N$$

  
From this system, we construct a sequence of values $\{u_{i}^{k}\}_{i=0,N}$ such that:

$$u_{i}^{k+1}=G_{i}(u_{0}^{k},u_{1}^{k},..u_{j}^{k}...,u_{N}^{k}), \forall i=1,N$$

  
If the sequence $\{u_{i}^{k}\}_{i=0,N}$ converges, it converges to a fixed point (i.e., a solution) and thus to the solution of the initial nonlinear problem. The convergence condition of the fixed-point sequence is given by the classical fixed-point theorem:

**Fixed-point theorem:**

The iterative sequence $(u^n)$ converges in the vicinity of a fixed point if and only if the Jacobian matrix $J$ of $G$: $J_{i,j}=\left(\frac{\partial G_{i}}{\partial u_{j}}\right)$, has a norm less than 1, i.e., has eigenvalues with magnitudes less than 1 in this vicinity.
  
### Linearized Implicit Scheme

$$\frac{u_{i}^{n+1}-u_{i}^{n}}{dt}-\left(\frac{\kappa_{i+\frac{1}{2}}^{n+1}\left(u_{i+1}^{n+1}-u_{i}^{n+1}\right)-\kappa_{i-\frac{1}{2}}^{n+1}\left(u_{i}^{n+1}-u_{i-1}^{n+1}\right)}{dx^{2}}\right)+\sigma\left(u_{i}^{n+1}\left(u_{i}^{n}\right)^{3}-1\right)=Q_{i}$$

This scheme can also be written as a fixed-point iteration:

$$\left[u_{i}^{n+1}\right]=\mathcal{A}^{-1}\left[u_{i}^{n}+dt (Q_{i}+\sigma)\right]$$

where $\mathcal{A}$ is the following tridiagonal matrix of order $N+1$:

$$\mathcal{A}=\left[\begin{array}{cccccc}
a_{0} & 2 b_0 & 0 & & 0 & 0\\
c_{1} & a_{1} & b_{1} & & 0 & 0\\
0 & c_{2} & a_{2} & \ddots & 0 & 0\\
& & \ddots & \ddots & \ddots & \\
0& 0 & 0 & \ddots & a_{N-1} & b_{N-1}\\
0 & 0 & 0 & & 0 & 1\end{array}\right]$$

with

$$a_{i}=1+dt\left(\frac{\kappa_{i+\frac{1}{2}}^{n}}{dx^{2}}+\frac{\kappa_{i-\frac{1}{2}}^{n}}{dx^{2}}+4\sigma\left(u_{i}^{n}\right)^{3}\right)$$
$$b_{i}=-dt\frac{\kappa_{i+\frac{1}{2}}^{n}}{dx^{2}}$$
$$c_{i}=-dt\frac{\kappa_{i-\frac{1}{2}}^{n}}{dx^{2}}$$

The convergence of the fixed-point iteration is related to the stability of the scheme. The Jacobian matrix $\mathcal{G}$ of the fixed-point iteration is written as:

$$\mathcal{G}_{i,j}=\left[\mathcal{A}_{i,j}+\frac{\partial\mathcal{A}_{i,j}}{\partial u_{k}}u_{k}^{n}\right]^{-1}$$

In the linear case, the matrix $\mathcal{A}$ is constant, and we have $\mathcal{G}=\mathcal{A}^{-1}$. The linear implicit scheme is unconditionally stable, and thus the matrix $\mathcal{A}^{-1}$ has eigenvalues with magnitudes less than 1. In the nonlinear case, the derivative of the coefficients of $\mathcal{A}$ with respect to the solution $u_{i}^{n}$ must be considered. The matrix $\frac{\partial\mathcal{A}}{\partial u_{k}}u_{k}^{n}$ being proportional to $dt$, for small time steps we have $\mathcal{G\approx A}^{-1}$, and the fixed-point iteration converges since the eigenvalues of $\mathcal{A}^{-1}$ are less than 1 in magnitude.

We can therefore conclude that the linearized implicit scheme converges for small time steps.

### Newton's Scheme

To find the roots of nonlinear equations, we can use the Newton-Raphson method, which consists of constructing the following fixed-point iterative sequence:

$$\displaystyle \left[u_{i}^{k+1}\right]=\left[u_{i}^{k}\right]-\left[\mathbf{J}^{k}\right]_{i,j}^{-1}\left[F_{j}(u_{0}^{k},..u_{N}^{k})\right]$$

$\mathbf{J}^{k}$ is the Jacobian matrix of the functions $F_{i}$: $\mathbf{J}_{i,j}^{k}=\frac{\partial F_{i}}{\partial u_{j}}(u_{0}^{k}..u_{N}^{k})$ calculated at iteration $k$. This relation is written in matrix form as:

$$\mathbf{J}^{k}\left[u_{i}^{k+1}-u_{i}^{k}\right]=-\left[F_{j}(u_{0}^{k},..u_{N}^{k})\right]$$

At each Newton iteration, a linear system of order $N+1$ must be solved.

In our case, the Jacobian matrix $\mathbf{J}^{k}$ is tridiagonal and is written as:

$$\displaystyle \mathbf{J}^{k}=\left[\begin{array}{cccccc}
a_{0} & 2 b_0 & 0 & & 0 & 0\\
c_{1} & a_{1} & b_{1} & & 0 & 0\\
0 & c_{2} & a_{2} & \ddots & 0 & 0\\
& & \ddots & \ddots & \ddots & \\
0& 0 & 0 & \ddots & a_{N-1} & b_{N-1}\\
0 & 0 & 0 & & 0 & 1\end{array}\right]$$

with:

$$a_{i} = \left(\frac{\kappa_{i+\frac{1}{2}}+\kappa_{i-\frac{1}{2}}}{dx^2}+\frac{\frac{1}{2}(\frac{\partial \kappa_i}{\partial u_i})(2u_i^k-u_{i-1}^k-u_{i+1}^k)}{dx^{2}}+4\sigma\left(u_{i}^k\right)^{3}\right),$$

$$b_{i} = -\frac{\kappa_{i+\frac{1}{2}}^k}{dx^{2}}+\frac{\frac{1}{2}(\frac{\partial\kappa_{i+1}}{\partial u_{i+1}})(u_{i+1}^k-u_{i}^k)}{dx^{2}}$$

$$c_{i} = -\frac{\kappa_{i-\frac{1}{2}}^k}{dx^{2}}+\frac{\frac{1}{2}(\frac{\partial\kappa_{i-1}}{\partial u_{i-1}})(u_{i-1}^k-u_{i}^k)}{dx^{2}}$$

As with the explicit scheme, we can simplify the calculation of these coefficients by neglecting the terms in $\frac{1}{2}(\frac{\partial\kappa_{i}}{\partial u_{i}})$ compared to those in $\kappa_{i}$, which gives:

$$a_{i}\approx\left(\frac{\kappa_{i+\frac{1}{2}}^{k}+\kappa_{i-\frac{1}{2}}^{k}}{dx^{2}}+4\sigma\left(u_{i}^{k}\right)^{3}\right)$$

$$b_{i}\approx-\frac{\kappa_{i+\frac{1}{2}}^{k}}{dx^{2}}$$

$$c_{i}\approx-\frac{\kappa_{i-\frac{1}{2}}^{k}}{dx^{2}}$$

# Implementation

Implement the linearized implicit scheme and the Newton method to solve the nonlinear diffusion problem.
We can use the following parameters:

1. $\kappa_{0}=0.01$, $\sigma=0.1$, $\beta = 1$, $\kappa(u)=\kappa_0 \sqrt u $,
2. $\kappa_{0}=0.01$, $\sigma=1$, $\beta = 300$, $\kappa(u)=\kappa_0 u^2$.

For the numerical application, we will take $N=50, ... 10000$ discretization points.

For the time step, we will take

$$dt=\gamma \frac{2}{4\sigma u_{max}^{3}+\frac{4\kappa(u_{max})}{dx^{2}}}$$

with $\gamma=0.1$, $1$ or $10$.

For the linear solver part, we will use the `hypre` library (https://github.com/hypre-space/hypre).
In particular, the IJ interface will be used to store the sparse matrix and the right-hand side vector (https://hypre.readthedocs.io/en/latest/ch-ij.html)

The use of a multigrid preconditioner can be considered by following the example provided in the `hypre` documentation (https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html).