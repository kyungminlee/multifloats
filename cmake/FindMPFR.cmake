# FindMPFR - locate GNU MPFR (libmpfr + mpfr.h) and its GMP dependency.
#
# Sets:
#   MPFR_FOUND
#
# Target:
#   MPFR::MPFR  (IMPORTED, includes gmp as a transitive link)

find_path(MPFR_INCLUDE_DIR mpfr.h
    HINTS ENV MPFR_DIR ENV MPFR_INCLUDE_DIR
    PATH_SUFFIXES include)

find_library(MPFR_LIBRARY NAMES mpfr
    HINTS ENV MPFR_DIR ENV MPFR_LIB_DIR
    PATH_SUFFIXES lib lib64)

find_library(GMP_LIBRARY NAMES gmp
    HINTS ENV GMP_DIR ENV GMP_LIB_DIR
    PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
    REQUIRED_VARS MPFR_LIBRARY MPFR_INCLUDE_DIR GMP_LIBRARY)

if(MPFR_FOUND AND NOT TARGET MPFR::MPFR)
    add_library(MPFR::MPFR UNKNOWN IMPORTED)
    set_target_properties(MPFR::MPFR PROPERTIES
        IMPORTED_LOCATION "${MPFR_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${GMP_LIBRARY}")
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY GMP_LIBRARY)
