project(OpenAT)
enable_language(Fortran)

set(OpenAT_SOURCES
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

set(OpenAT_OBJECTS "")
foreach(src IN LISTS OpenAT_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}_obj)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND OpenAT_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

add_library(OpenAT STATIC ${OpenAT_OBJECTS})