cmake_minimum_required(VERSION 3.5)
project(bicgstab)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_C_COMPILER "/usr/local/bin/icc" CACHE STRING "C compiler" FORCE)
set(CMAKE_C_FLAGS "-O3 -m64 -fopenmp -shared-intel -lifport -lifcore" CACHE STRING "C flags" FORCE)
set(CMAKE_C_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "C linker flags" FORCE)
set(CMAKE_FORTRAN_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "CXX linker flags" FORCE)

message(STATUS "CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")

find_package(XabClib REQUIRED)

if(XabClib_FOUND)
    message(STATUS "XabClib found")
else()
    message(FATAL_ERROR "XabClib not found")
endif()

add_executable(${PROJECT_NAME}
    sampleBICGSTAB.c
)

# link fortran static library
target_link_libraries(${PROJECT_NAME}
    PRIVATE XabClib::XabClib
)
