cmake_minimum_required(VERSION 3.6)
project(ananth_sampling)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp BoundedTableaux.cpp Tableaux.cpp BoundedTableaux.h Tableaux.h SimpleTableaux.cpp SimpleTableaux.h LozengeTiling.cpp LozengeTiling.h Coupling.cpp Coupling.h)
add_executable(ananth_sampling ${SOURCE_FILES})