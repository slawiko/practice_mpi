cmake_minimum_required(VERSION 3.8)
project(MPI)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        main.c)

add_executable(MPI ${SOURCE_FILES})