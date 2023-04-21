project(LANCZOS)
enable_language(Fortran)

set(LANCZOS_SOURCES
    BLAS_dstevd.f90
    LAPACK_dstevd.f90
    Xabclib_LANCZOS.f90
    Xabclib_MGramSch.f90
    Xabclib_QSORTD.f90
)

set(LANCZOS_OBJECTS "")
foreach(src IN LISTS LANCZOS_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}_obj)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND LANCZOS_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

# set(LANCZOS_OBJECTS ${LANCZOS_OBJECTS} PARENT_SCOPE)
add_library(LANCZOS STATIC ${LANCZOS_OBJECTS})