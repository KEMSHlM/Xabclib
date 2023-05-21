cmake_minimum_required(VERSION 3.5)

project(XabClib VERSION 1.0.0)
enable_language(Fortran)

# intel compiler
if(APPLE)
    # icc compiler (m2 mac)
    set(CMAKE_Fortran_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler flags" FORCE)
    set(CMAKE_Fortran_FLAGS "-O3 -m64 -fopenmp -shared-intel -lifport -lifcore -march=native -diag-disable=10441" CACHE STRING "Fortran compiler flags" FORCE)
elseif(UNIX)
    # icc compiler (ito)
    set(CMAKE_Fortran_COMPILER "/home/app/intel/oneapi/2021.3/compiler/2021.3.0/linux/bin/intel64/ifort" CACHE STRING "Fortran compiler flags" FORCE)
    set(CMAKE_Fortran_FLAGS "-O3 -m64 -fopenmp -shared-intel -lifport -lifcore -march=native -diag-disable=10441" CACHE STRING "Fortran compiler flags" FORCE)
endif()

message(STATUS "CMAKE_Fortran_COMPILER: ${CMAKE_FORTRAN_COMPILER}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_FORTRAN_FLAGS}")

# サブディレクトリ内のCMakeLists.txtファイルを追加
add_subdirectory(OpenAT)
add_subdirectory(Xabclib/Xabclib_GMRES)
add_subdirectory(Xabclib/Xabclib_BICGSTAB)
add_subdirectory(Xabclib/Xabclib_LANCZOS)
add_subdirectory(Xabclib/Xabclib_ARNOLDI)
add_subdirectory(Xabclib/Xabclib_CG)

# Create the XabClib library
add_library(XabClib)

# オブジェクトファイルをまとめて静的ライブラリにする
target_link_libraries(XabClib PUBLIC
    OpenAT
    GMRES
    BICGSTAB
    LANCZOS
    ARNOLDI
    CG
)

# Specify the install location
install(TARGETS ${PROJECT_NAME}
    EXPORT XabClibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)

# Create the MyLibTargets.cmake in the build directory
export(EXPORT XabClibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/XabClib/XabClibTargets.cmake"
    NAMESPACE XabClib::
)

# Install the export set for use with the install-tree
install(EXPORT XabClibTargets
    FILE XabClibTargets.cmake
    NAMESPACE XabClib::
    DESTINATION lib/cmake/XabClib
)

# Create a Config.cmake file for find_package()
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/XabClib/XabClibConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/CMake/XabClibConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/XabClib/XabClibConfig.cmake"
    INSTALL_DESTINATION lib/cmake/XabClib
)

# Install the Config files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/XabClib/XabClibConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/XabClib/XabClibConfigVersion.cmake"
        DESTINATION lib/cmake/XabClib)