project(BICGSTAB)
enable_language(Fortran)

set(BICGSTAB_SOURCES
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

set(BICGSTAB_OBJECTS "")
foreach(src IN LISTS BICGSTAB_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}_obj)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND BICGSTAB_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

# set(BICGSTAB_OBJECTS ${BICGSTAB_OBJECTS} PARENT_SCOPE)
add_library(BICGSTAB STATIC ${BICGSTAB_OBJECTS})