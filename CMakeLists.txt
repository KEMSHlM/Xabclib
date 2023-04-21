cmake_minimum_required(VERSION 3.5)

project(XabClib)
enable_language(Fortran)

# intel compiler
set(CMAKE_FORTRAN_COMPILER "/opt/intel/oneapi/compiler/2023.0.0/mac/bin/intel64/ifort" CACHE STRING "Fortran compiler" FORCE)
set(CMAKE_FORTRAN_FLAGS "-O0 -m64 -fopenmp -shared-intel -cpp -fixed" CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_FORTRAN_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "FORTRAN linker flags" FORCE)
set(CMAKE_C_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "C linker flags" FORCE)

# サブディレクトリ内のCMakeLists.txtファイルを追加
add_subdirectory(OpenAT)
add_subdirectory(Xabclib/Xabclib_GMRES)
add_subdirectory(Xabclib/Xabclib_BICGSTAB)
add_subdirectory(Xabclib/Xabclib_LANCZOS)
add_subdirectory(Xabclib/Xabclib_ARNOLDI)
add_subdirectory(Xabclib/Xabclib_CG)

# オブジェクトファイルをまとめて静的ライブラリにする
file(WRITE dummy.c "")
add_library(XabClib STATIC dummy.c)

target_link_libraries(XabClib PUBLIC
    OpenAT
    GMRES
    BICGSTAB
    LANCZOS
    ARNOLDI
    CG
)