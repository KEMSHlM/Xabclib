project(LANCZOS)
enable_language(Fortran)

add_library(LANCZOS OBJECT
    BLAS_dstevd.f90
    LAPACK_dstevd.f90
    Xabclib_LANCZOS.f90
    Xabclib_MGramSch.f90
    Xabclib_QSORTD.f90
)

install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
