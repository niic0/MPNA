#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    {
      for(int j = 0; j<m; j++) {
	A[i][j] = At[j][i];
      }
    }
}

void fill_matrix(int n, int m, double A[n][m])
{
  for (int i = 0; i < n; i++)
    {
      for (int k = 0; k < m; k++)
        {
	  A[i][k] = 1.0;
        }
    }
}


int main(int argc, char *argv[])
{  
  // 100 time measures
  int number_of_measure = 10;
  double time_taken[number_of_measure];
  int N = 100;
  
  printf("Working with square matrices with variable size\n");

  // Open file to save data for ploting
  FILE *data_file = fopen("data.txt", "w");
  if (!data_file) {
    perror("Error creating data file");
    return 1;
  }
  fprintf(data_file, "# N Time_Taken\n");

  for(int i = 0; i<number_of_measure; i++) {
    // Allocate memory for matrices
    // Note how we allocate memory for a 2D array in a single block
    // This is done to ensure that the memory is contiguous
    // Can you comment on how data are stored in memory for this 2D array?
    double **A = (double **)malloc(N * N * sizeof(double));
    double **B = (double **)malloc(N * N * sizeof(double));
    double **Bt = (double **)malloc(N * N * sizeof(double)); // B transpose
    double **C = (double **)malloc(N * N * sizeof(double));

    // Initialize matrices
    fill_matrix(N, N, A);
    fill_matrix(N, N, B);
    fill_matrix(N, N, C);

    transpose_matrix(N, N, B, Bt);
    
    // Call matrix_matrix_multiplication ijk
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    matrix_matrix_multiplication(N, N, N, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      
    time_taken[i] = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("[ijk] N=%d ; Time taken: %e s\n", N, time_taken[i]);
      
    // Write data to file
    fprintf(data_file, "%d %f\n", N, time_taken);

    N += 100;
    // Free memory
    free(A);
    free(B);
    free(C);
  }

  return EXIT_SUCCESS;
}
