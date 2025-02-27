cmake_minimum_required(VERSION 3.13)
project(SparseSolvers)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Optimization options
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -funroll-loops")

# Include path
include_directories(include)

# Sources files
set(SOURCES
    src/ConjugateGradient.cpp
    src/GaussSeidel.cpp
    src/Jacobi.cpp
    src/Utils.cpp
    src/Benchmark.cpp
    src/ArgumentParser.cpp
    src/GMRES.cpp
    main.cpp
)

# Executable creation
add_executable(SparseSolvers ${SOURCES})