cmake_minimum_required(VERSION 3.19)
project(cpp_linalg)

set(CMAKE_CXX_STANDARD 14)

#add_library(include/cpp_linalg/mem_blocks.h include/cpp_linalg/matrix_core.h)
add_executable(cpp_linalg src/app.cpp )#include/cpp_linalg/mem_blocks.h include/cpp_linalg/matrix_core.h )