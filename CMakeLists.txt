cmake_minimum_required(VERSION 3.5)

project(onigiri)
enable_language(Fortran)
set(CMAKE_FORTRAN_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler" FORCE)
set(CMAKE_FORTRAN_FLAGS "-O3 -O3 -m64 -fopenmp -shared-intel -cpp -fixed" CACHE STRING "Fortran compiler flags" FORCE)

add_subdirectory(XabClib)
add_subdirectory(OpenAT)

add_library(${PROJECT_NAME}
    ${XABCLIB_SRCS}
    ${OPENAT_SRCS}
)
