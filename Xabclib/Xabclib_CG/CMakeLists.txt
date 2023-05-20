project(CG)
enable_language(Fortran)

add_library(CG OBJECT
    Xabclib_CG.f90
    Xabclib_CG_PCHK.f90
    Xabclib_DSRMV_ForJacobiIter.f90
    Xabclib_EVAL_SYM.f90
    Xabclib_IC0D.f90
    Xabclib_PCGENS.f90
    Xabclib_PCSLV_SYM.f90
)

install(TARGETS ${PROJECT_NAME}
    EXPORT XabClibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
