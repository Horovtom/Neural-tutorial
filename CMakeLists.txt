cmake_minimum_required(VERSION 3.1)
project(Neural-net-tutorial)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES Neural-net-tutorial/Neural-net-tutorial.cpp)
add_executable(NeuralNet ${SOURCE_FILES})