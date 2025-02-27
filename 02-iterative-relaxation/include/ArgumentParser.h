#pragma once

#include <string>

/**
 * @brief Stores user-defined parameters for solver execution.
 */
struct SolverConfig {
    std::string input_file = "";
    size_t max_iterations = 100000;
    double tolerance = 1e-6;
    int mesh_size = -1;  ///< If set, use Poisson matrix instead of file.
};

void parseArguments(int argc, char *argv[], SolverConfig &config);