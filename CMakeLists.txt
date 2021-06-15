cmake_minimum_required(VERSION 3.19)
project(cpp_linalg)

set(CMAKE_CXX_STANDARD 14)

add_library(include/cpp_linalg/matrix.h include/cpp_linalg/mem_blocks.h)
#add_executable(cpp_linalg_test_app src/app.cpp src/core.cpp include/cpp_linalg/core.h src/matrix.cpp src/mem_blocks.cpp)