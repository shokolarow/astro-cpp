cmake_minimum_required(VERSION 3.17)
project(astro_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O3 -mfma -mavx2")

include_directories(src)
include_directories(/home/alexshokolarov/Documents/Resources/cpplibs/)

add_executable(astro_cpp
        src/constants.h
        src/main.cpp
        src/Propagator.h src/Propagator.cpp src/tools.cpp src/tools.h)
