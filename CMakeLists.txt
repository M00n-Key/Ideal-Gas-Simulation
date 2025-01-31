
cmake_minimum_required(VERSION 3.24)

project(IdealGasSimulation)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${IdealGasSimulation_SOURCE_DIR}/bin)

include_directories(${IdealGasSimulation_SOURCE_DIR}/include)

add_executable(simulation src/main.cpp)


