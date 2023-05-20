project(OpenAT)
enable_language(Fortran)

add_library(OpenAT OBJECT
    OpenATI_DAFGS.f90
    OpenATI_DAFRT.f90
    OpenATI_DSRMV.f90
    OpenATI_DSRMV_Setup.f90
    OpenATI_DURMV.f90
    OpenATI_DURMV_Setup.f90
    OpenATI_EIGENSOLVE.f90
    OpenATI_LINEARSOLVE.f90
    OpenATI_NZDIV.f90
    OpenATI_UTIL.f90
    OpenATI_init.f90
    OpenATI_GET_VER.f90
    OpenATI_DAFSTG.f90
    OpenATI_DAFMC_CCS2CRS.f90
    OpenATI_QSORTI.f90
)

install(TARGETS ${PROJECT_NAME}
    EXPORT XabClibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
