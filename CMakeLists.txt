cmake_minimum_required(VERSION 3.5)

project(XabClib VERSION 1.0.0)
enable_language(Fortran)

# intel compiler
find_program(IFORT_COMPILER ifort)
if(IFORT_COMPILER)
    message(STATUS "ifort compiler found at ${IFORT_COMPILER}")
    set(CMAKE_Fortran_COMPILER ${IFORT_COMPILER} CACHE STRING "Fortran compiler" FORCE)
    set(CMAKE_Fortran_FLAGS "-O3 -m64 -shared-intel -lifport -lifcore -march=native -fixed" CACHE STRING "Fortran compiler flags" FORCE)
else()
    message(FATAL_ERROR "ifort compiler not found")
endif()

message("-- Searching OpenMP Library")
# intel OpenMP library path
if (APPLE)
    set(OpenMP_PATH "$ENV{ONEAPI_ROOT}/compiler/latest/mac/compiler/lib")
elseif (UNIX)
    set(OpenMP_PATH "$ENV{ONEAPI_ROOT}/compiler/latest/linux/compiler/lib/intel64")
elseif (WIN32)
    # 僕は使わないから設定してない
    # set(OpenMP_PATH "C:/path/for/windows")
else()
    message(FATAL_ERROR "Unsupported operating system")
endif()

find_library(OPENMP_LIB 
    NAMES iomp5
    HINTS ${OpenMP_PATH}
)

if (OPENMP_LIB)
    message(STATUS "OpenMP found : ${OPENMP_LIB}")    
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
    message(STATUS "CMAKE_Fortran_FLAGS : ${CMAKE_Fortran_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found")    
endif()

message(STATUS "CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")

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