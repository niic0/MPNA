# Matrix-Matrix Multiplication

## Description

The goal of this exercise is to implement a function that multiplies two matrices.

Given two matrices `A` and `B` of dimensions `m x n` and `n x p` respectively, the product of `A` and `B` is a matrix
`C` of dimensions `m x p` such that `C[i][j]` is the dot product of the `i`-th row of `A` and the `j`-th column of `B`.

In the following, we will implement in C such a function. The matrices will be represented as two-dimensional arrays of
doubles.

To compile using CMake:

```bash
cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug .
cmake --build build-debug
```

or, in release mode:

```bash
cmake -B build-release -DCMAKE_BUILD_TYPE=Release .
cmake --build build-release
```

You can then run the program with:

```bash
./build-debug/matrix_matrix <n> <m> <p>
```

where `n`, `m`, and `p` are the dimensions of the matrices.

## Exercises

### Exercise 1

Write one function `matrix_matrix_multiplication` that computes a matrix-matrix products. It should take as input these
arguments:

- the dimensions `m`, `n`, and `p` of the matrices.
- three matrices `A`, `B`, and `C` of dimensions `m x n`, `n x p`, and `m x p` respectively.

Is the function working correctly? How can you test it?

### Exercise 2

Compute the complexity of the function `matrix_matrix_multiplication`.

Can this function written differently? If so, how?

### Exercise 3

We want to benchmark the function `matrix_matrix_multiplication`. For the case of square matrices `m = n = p`, generate
two random matrices `A` and `B` and compute their product. Use the `time` command to measure the time it takes to
compute the product.

Repeat the experiment for different values of `n` and plot the time it takes to compute the product as a function of
`n`.

Compare the results with those of the alternative implementations you have proposed in Exercise 2.

### Exercise 4

Matrices data structures are key in many scientific applications. In the case of dense matrices, the two-dimensional
array representation is the most common.

Provide a specialized implementation of `matrix_matrix_multiplication` for the case of square matrices which operates on
`A`, `B transposed`, and `C`. Compare the performance of this implementation with the generic one.

### Exercise 5
> Compare the performance of your implementations with the one provided by the BLAS library. You can use the `cblas_dgemm`
or `dgemm` functions from the BLAS library.
