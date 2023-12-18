cmake_minimum_required(VERSION 3.0)

project(sim)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/json/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/json/single_include")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -L/usr/lib/x86_64-linux-gnu")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -L/usr/lib/x86_64-linux-gnu")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(run src/main.cpp src/sim.cpp)