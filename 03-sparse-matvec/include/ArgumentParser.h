#pragma once

#include <string>

/**
 * @brief Stores user-defined parameters for solver execution.
 */
struct CliConfig {
    std::string input_file = "";
    size_t mpi_node = 1;
};

void parseArguments(int argc, char *argv[], CliConfig &config);