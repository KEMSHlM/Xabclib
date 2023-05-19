project(ARNOLDI)
enable_language(Fortran)

add_library(ARNOLDI OBJECT
    BLAS_dgeev.f90
    LAPACK_dgeev.f90
    Xabclib_ARNOLDI.f90
)

install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
