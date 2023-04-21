cmake_minimum_required(VERSION 3.5)

project(XabClib)
enable_language(Fortran)

# intel compiler
set(CMAKE_FORTRAN_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler" FORCE)
set(CMAKE_FORTRAN_FLAGS "-O0 -m64 -fopenmp -shared-intel -cpp -fixed" CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_FORTRAN_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "FORTRAN linker flags" FORCE)
set(CMAKE_C_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "C linker flags" FORCE)

# gcc-12 compiler
# set(CMAKE_FORTRAN_COMPILER "/opt/homebrew/bin/gfortran" CACHE STRING "Fortran compiler" FORCE)
# set(CMAKE_FORTRAN_FLAGS "-Os -m64 -fopenmp -cpp -fixed -DHF_COMPILER -march=armv8-a" CACHE STRING "Fortran compiler flags" FORCE)
message(STATUS "CMAKE_FORTRAN_COMPILER = ${CMAKE_FORTRAN_COMPILER}")

add_subdirectory(Xabclib)
add_subdirectory(OpenAT)

add_library(${PROJECT_NAME} STATIC
    ${XABCLIB_SRCS}
    ${OPENAT_SRCS}
)