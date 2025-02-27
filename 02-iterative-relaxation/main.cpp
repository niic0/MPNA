#include "include/ArgumentParser.h"
#include "include/Benchmark.h"
#include "include/CSRMatrix.h"
#include "include/Solver.h"
#include "include/Utils.h"

int main(int argc, char *argv[]) {
    SolverConfig config;
    parseArguments(argc, argv, config);  // Handle user input

    // Load matrix
    CSRMatrix A;
    if (!config.input_file.empty()) {
        A.fillWithFile(config.input_file);
    } else if (config.mesh_size > 0) {
        A.fillWithPoisson(config.mesh_size);
    } else {
        std::cerr << "Error: No input file or mesh size specified. Use -h for help.\n";
        return 1;
    }

    // Initialize vectors
    std::vector<double> b(A.getRows(), 1.0);
    std::vector<double> x(A.getRows(), 0.0);
    std::vector<SolverResult> results;

    // Run solvers
    results.push_back(runSolver("Jacobi", jacobi, A, b, x,
                                config.max_iterations, config.tolerance));
    results.push_back(runSolver("Gauss-Seidel", gaussSeidel, A, b, x,
                                config.max_iterations, config.tolerance));
    results.push_back(runSolver("Conjugate Gradient", conjugateGradient, A, b,
                                x, config.max_iterations, config.tolerance));
    results.push_back(runSolver("GMRES", gmres, A, b,
                                x, config.max_iterations, config.tolerance));

    // Print results
    printResultsTable(results);

    return 0;
}
