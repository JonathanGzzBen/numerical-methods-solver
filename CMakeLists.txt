cmake_minimum_required(VERSION 3.14)

project(numerical-methods-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src) # Library

add_subdirectory(console)
