cmake_minimum_required(VERSION 3.2)
project(quadtree)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -fopenmp -ffast-math")

set(SOURCE_FILES main.cpp quadtree.cpp quadtree.h treecode.cpp treecode.h)

add_executable(quadtree ${SOURCE_FILES})
target_link_libraries(quadtree m)
