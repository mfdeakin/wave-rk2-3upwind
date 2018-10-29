
cmake_minimum_required (VERSION 3.3)
project(waves)

include_directories(./include)

add_library(waves ./src/wave_eqn_small.cpp ./src/wave_eqn_large.cpp)
add_executable(waves_tester ./src/test.cpp)
add_executable(waves_solver ./src/solver.cpp)
add_test(waves_test waves_tester)

target_link_libraries(waves_tester gtest waves)
target_link_libraries(waves_solver waves)
