cmake_minimum_required(VERSION 3.13)
project(NonlinearDiffusionSolver)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Optimization Flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -funroll-loops")

# Enable Ninja
set(CMAKE_GENERATOR "Ninja")

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Find HYPRE (Fix Missing Library Issue)
set(HYPRE_DIR ${CMAKE_SOURCE_DIR}/external/install)
include_directories(${HYPRE_DIR}/include)
link_directories(${HYPRE_DIR}/lib)

# Sources
set(SOURCES
    src/HypreInterface.cpp
    src/Utils.cpp
    src/Solver.cpp
    src/MatrixBuilder.cpp
    src/main.cpp
)

# Create Executable
add_executable(diffusion ${SOURCES})

# Correctly Link HYPRE
target_link_libraries(diffusion PRIVATE ${HYPRE_DIR}/lib/libHYPRE.a ${MPI_CXX_LIBRARIES})
