cmake_minimum_required(VERSION 3.5)
project(bicgstab LANGUAGES C)

set(CMAKE_C_COMPILER "/usr/local/bin/icc" CACHE STRING "C compiler" FORCE)
set(CMAKE_C_FLAGS "-O3 -m64 -fopenmp -shared-intel -lifport -lifcore" CACHE STRING "C flags" FORCE)

message(STATUS "CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")

find_package(XabClib 1.0.0 REQUIRED)

if(XabClib_FOUND)
    message(STATUS "XabClib is found")
else()
    message(FATAL_ERROR "XabClib is not found")
endif()

add_executable(${PROJECT_NAME}
    sampleBICGSTAB.c
)

# link fortran static library
target_link_libraries(${PROJECT_NAME}
    PRIVATE XabClib::XabClib
)
