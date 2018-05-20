cmake_minimum_required(VERSION 3.6)
project(ananth_sampling)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Tableaux.cpp Tableaux.h)
add_executable(ananth_sampling ${SOURCE_FILES})