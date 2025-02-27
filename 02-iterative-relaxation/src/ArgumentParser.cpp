#include "../include/ArgumentParser.h"
#include <getopt.h> // POSIX-style argument parser
#include <iostream>

/**
 * @brief Displays help message for command-line usage.
 */
void printHelp() {
  std::cout
      << "Usage: ./SparseSolver [OPTIONS]\n"
      << "Options:\n"
      << "  -i, --input <file>       Path to input matrix file (.mtx)\n"
      << "  -m, --mesh <size>        Generate a Poisson matrix with given mesh "
         "size\n"
      << "  -n, --iterations <num>   Set maximum iterations (default: 100000)\n"
      << "  -t, --tolerance <value>  Set convergence tolerance (default: "
         "1e-6)\n"
      << "  -h, --help               Display this help message\n";
}

/**
 * @brief Parses command-line arguments and fills SolverConfig.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param config Reference to SolverConfig to store parsed values.
 */
void parseArguments(int argc, char *argv[], SolverConfig &config) {
  static struct option long_options[] = {
      {"input", required_argument, nullptr, 'i'},
      {"mesh", required_argument, nullptr, 'm'},
      {"iterations", required_argument, nullptr, 'n'},
      {"tolerance", required_argument, nullptr, 't'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "i:m:n:t:h", long_options, nullptr)) !=
         -1) {
    switch (opt) {
    case 'i':
      config.input_file = optarg;
      break;
    case 'm':
      config.mesh_size = std::stoi(optarg);
      break;
    case 'n':
      config.max_iterations = std::stoi(optarg);
      break;
    case 't':
      config.tolerance = std::stod(optarg);
      break;
    case 'h':
      printHelp();
      exit(0);
    default:
      std::cerr << "Invalid argument. Use -h or --help for usage.\n";
      exit(1);
    }
  }
}