cmake_minimum_required(VERSION 3.5)

project(XabClib)
enable_language(Fortran)
set(CMAKE_FORTRAN_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler" FORCE)
set(CMAKE_FORTRAN_FLAGS "-O3 -O3 -m64 -fopenmp -shared-intel -cpp -fixed" CACHE STRING "Fortran compiler flags" FORCE)

add_subdirectory(XabClib)
add_subdirectory(OpenAT)

add_library(${PROJECT_NAME}
    ${XABCLIB_SRCS}
    ${OPENAT_SRCS}
)

# インストールの設定
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

export(TARGETS ${PROJECT_NAME}
       NAMESPACE ${PROJECT_NAME}::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
)

include(CMakePackageConfigHelpers)

write_basic_config_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION 1.0.0
    COMPATIBILITY AnyNewerVersion
)

install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION lib/cmake/${PROJECT_NAME}
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME}
)
