cmake_minimum_required(VERSION 3.5)

project(XabClib LANGUAGES Fortran)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# intel compiler
set(CMAKE_FORTRAN_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler" FORCE)
set(CMAKE_FORTRAN_FLAGS "-O3 -m64 -fopenmp -shared-intel -cpp -fixed" CACHE STRING "Fortran compiler flags" FORCE)
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

# インストールの設定
# プロジェクトのインストール先への内容を指定する
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# 上のターゲット設定XabClibTargets.cmakeという名前で出力する
export(TARGETS ${PROJECT_NAME}
       NAMESPACE ${PROJECT_NAME}::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
)

# パッケージ設定ファイルを作成するための関数を提供する
include(CMakePackageConfigHelpers)

# パッケージ設定ファイルを出力する
write_basic_config_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION 1.0.0
    COMPATIBILITY AnyNewerVersion
)

# 出力されたバージョン設定ファイルを適切なディレクトリにインストールする
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION lib/cmake/${PROJECT_NAME}
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME}
)