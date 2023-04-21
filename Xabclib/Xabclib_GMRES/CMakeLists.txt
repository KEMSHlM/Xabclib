project(GMRES)
enable_language(Fortran)

set(GMRES_SOURCES
    Xabclib_GMRES.f90
    Xabclib_GMRES2.f90
    Xabclib_GMRES_GS_Policy.f90
    Xabclib_GMRES_INIT.f90
    Xabclib_GMRES_MGS_MIN.f90
    Xabclib_GMRES_MIN.f90
    Xabclib_GMRES_PCHK.f90
)

set(GMRES_OBJECTS "")
foreach(src IN LISTS GMRES_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}_obj)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND GMRES_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

# set(GMRES_OBJECTS ${GMRES_OBJECTS} PARENT_SCOPE)
add_library(GMRES STATIC ${GMRES_OBJECTS})