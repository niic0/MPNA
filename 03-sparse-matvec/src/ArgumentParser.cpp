#include "../include/ArgumentParser.h"
#include <getopt.h> // POSIX-style argument parser
#include <iostream>

/**
 * @brief Displays help message for command-line usage.
 */
void printHelp() {
  std::cout << "Usage: ./matvec [OPTIONS]\n"
            << "Options:\n"
            << "  -i, --input <file>       Path to input matrix file (.mtx)\n"
               "size\n"
            << "  -h, --help               Display this help message\n";
}

/**
 * @brief Parses command-line arguments and fills SolverConfig.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param config Reference to CliConfig to store parsed values.
 */
void parseArguments(int argc, char *argv[], CliConfig &config) {
  static struct option long_options[] = {
      {"input", required_argument, nullptr, 'i'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "i:h", long_options, nullptr)) !=
         -1) {
    switch (opt) {
    case 'i':
      config.input_file = optarg;
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