project(BICGSTAB)
enable_language(Fortran)

add_library(BICGSTAB OBJECT
    ILUT.f90
    ILUTSOL.f90
    Xabclib_BICGSTAB.f90
    Xabclib_BICGSTABV.f90
    Xabclib_BICGSTAB_PCHK.f90
    Xabclib_BLAS1.f90
    Xabclib_EVAL.f90
    Xabclib_IDIAGPGEN.f90
    Xabclib_PCGENE.f90
    Xabclib_PCSLV.f90
)

install(TARGETS ${PROJECT_NAME}
    EXPORT XabClibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
