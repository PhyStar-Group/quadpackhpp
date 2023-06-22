
cmake_minimum_required (VERSION 3.8)

option(USE_CLANG "use latest clang compiler" OFF)
if (USE_CLANG)
    set(CMAKE_CXX_COMPILER clang++-10)
    set(CMAKE_C_COMPILER clang-10)
else()
    set(CMAKE_CXX_COMPILER g++-10)
    set(CMAKE_C_COMPILER gcc-10)
endif()

set(CMAKE_CXX_STANDARD 20)
project ("quadpack")


aux_source_directory(. DIR_SRCS)
aux_source_directory(./test DIR_SRCS)
add_executable (quadpack ${DIR_SRCS})

