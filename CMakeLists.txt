
cmake_minimum_required (VERSION 3.3)
project(waves)

enable_testing()

include_directories(./include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(pybind11 CONFIG REQUIRED)
include_directories(${pybind11_INCLUDE_DIRS})

add_library(waves ./src/wave_eqn_small.cpp ./src/wave_eqn_med.cpp ./src/wave_eqn_large.cpp)
add_executable(waves_tester ./src/test.cpp)
add_executable(waves_solver ./src/solver.cpp)
add_test(waves_test waves_tester)

target_link_libraries(waves_tester gtest waves)
target_link_libraries(waves_solver waves pybind11::embed)
