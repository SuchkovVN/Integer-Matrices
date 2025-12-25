cmake_minimum_required(VERSION 3.21)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_COMPILER "clang-20")
set(CMAKE_CXX_COMPILER "clang++-20")

project(IntegerMatrices VERSION 1.0.0 LANGUAGES C CXX)

find_library(PARI_LIB pari REQUIRED)

add_subdirectory(pari_wrapper)
add_subdirectory(utils)
add_subdirectory(compute)
