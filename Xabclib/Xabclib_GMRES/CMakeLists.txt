project(GMRES)
enable_language(Fortran)

add_library(GMRES OBJECT
    Xabclib_GMRES.f90
    Xabclib_GMRES2.f90
    Xabclib_GMRES_GS_Policy.f90
    Xabclib_GMRES_INIT.f90
    Xabclib_GMRES_MGS_MIN.f90
    Xabclib_GMRES_MIN.f90
    Xabclib_GMRES_PCHK.f90
)

install(TARGETS ${PROJECT_NAME}
    EXPORT XabClibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)