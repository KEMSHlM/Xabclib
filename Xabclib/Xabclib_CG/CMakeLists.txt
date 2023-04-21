project(CG)
enable_language(Fortran)

set(CG_SOURCES
    Xabclib_CG.f90
    Xabclib_CG_PCHK.f90
    Xabclib_DSRMV_ForJacobiIter.f90
    Xabclib_EVAL_SYM.f90
    Xabclib_IC0D.f90
    Xabclib_PCGENS.f90
    Xabclib_PCSLV_SYM.f90
)

set(CG_OBJECTS "")
foreach(src IN LISTS CG_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}.o)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND CG_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

# set(CG_OBJECTS ${CG_OBJECTS} PARENT_SCOPE)
add_library(CG STATIC ${CG_OBJECTS})