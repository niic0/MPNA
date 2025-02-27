cmake_minimum_required(VERSION 3.13)
project(matvec)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -funroll-loops")

include_directories(include)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(SOURCES
    src/Utils.cpp
    src/ArgumentParser.cpp
    src/MatrixVector.cpp
    main.cpp
)

add_executable(matvec ${SOURCES})
target_link_libraries(matvec MPI::MPI_CXX)
