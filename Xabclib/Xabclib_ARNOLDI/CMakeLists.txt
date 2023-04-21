project(ARNOLDI)
enable_language(Fortran)

set(ARNOLDI_SOURCES
    BLAS_dgeev.f90
    LAPACK_dgeev.f90
    Xabclib_ARNOLDI.f90
)

set(ARNOLDI_OBJECTS "")
foreach(src IN LISTS ARNOLDI_SOURCES)
    get_filename_component(name ${src} NAME_WE)
    set(name_obj ${name}_obj)
    add_library(${name_obj} OBJECT ${src})
    list(APPEND ARNOLDI_OBJECTS $<TARGET_OBJECTS:${name_obj}>)
endforeach()

# set(ARNOLDI_OBJECTS ${ARNOLDI_OBJECTS} PARENT_SCOPE)
add_library(ARNOLDI STATIC ${ARNOLDI_OBJECTS})