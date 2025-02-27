# Exerise 1
```c
void matrix_matrix_multiplication(int m, int n, int p, double A[m][n], double B[n][p], double C[m][p])
{
  for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
		  C[i][j] = 0.0;
		  for (int k = 0; k < p; k++) {
			  C[i][j] += A[i][k] * A[k][j];
	      }
	   }
   }
}
```

Compile with flag `-Wno-incompatible-pointer-types` flag with gcc: `gcc -Wno-incompatible-pointer-types main.c`.

> Is the function working correctly? How can you test it?

The function can be tested with different method:
- The result of the matrix multiply function can be compared to the result of a well-known library such as blas for example.
- The result can also be compared to a matrix that we the already calculated. This method implies that the initial matrix is pass as an argument.

# Exercise 2
> Compute the complexity of the function `matrix_matrix_multiplication`.
The complexity of GEMM is O(N^3) in this case.

> Can this function written differently? If so, how?
Yes, we can optimize this function in many ways. The performance can be enhance with different compiler. Compilers tend to know we are doing a matrix multiply and directly optimize is. But we can optimize the order in which arrays are loaded, by permuting indexes i, j and k.

Different time with different indexes order:
```
ijk: Time taken: 1.638484e-02 s
ikj: Time taken: 1.438192e-02 s
kij: Time taken: 1.437186e-02 s
kji: Time taken: 1.638166e-02 s
jik: Time taken: 1.672594e-02 s
jki: Time taken: 1.647758e-02 s
```

# Exercise 3
> Repeat the experiment for different values of `n` and plot the time it takes to compute the product as a function of `n`.

The program uses GNU plot and generates a data.txt. You can generate the image with the comma,nd:

```bash
$ gnuplot plot_script.gnu
```

All images with different optimization are stored in `plot` folder.

# Exercise 4
> Provide a specialized implementation of `matrix_matrix_multiplication` for the case of square matrices which operates on `A`, `B` transposed, and `C`. Compare the performance of this implementation with the generic one.

```c
void matrix_matrix_multiplication(int m, int n, int p, double A[m][n], double B[n][p], double C[m][p])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
		  C[i][j] += A[i][k] * B[j][k];
}

void transpose_matrix(int m, int n, double A[m][n], double At[n][m])
{
  for(int i = 0; i<m; i++)
      for(int j = 0; j<m; j++)
		  A[i][j] = At[j][i];
}
```

# Exercise 5
> Compare the performance of your implementations with the one provided by the BLAS library. You can use the cblas_dgemm or dgemm functions from the BLAS library.

```c
#include <gsl/gsl_cblas.h>
...
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, A, N, B, N, 1, C, N);
```

Compilation:
```bash
$ gcc -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm  -Wno-incompatible-pointer-types main.c
```

Résultat dans `performances_plot.png`


## Ressources
`BLIS`: Implémentation de BLAS pour les CPUs AMD. Des vidéos expliques comment à partir d'un problème général donné, arrivé ç construire une solution qui utlise des briques de BLIS (DGEMM, et autre problème d'algèbre linéaire).
