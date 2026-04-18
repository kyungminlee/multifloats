project(multifloats LANGUAGES CXX)

include(GNUInstallDirs)

add_library(multifloats STATIC multifloats_math.cc)
target_include_directories(multifloats
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# ABI bookkeeping — kept in sync with MULTIFLOATS_ABI_VERSION in multifloats_c.h.
# SOVERSION is only honored on shared libraries; it is set anyway so a future
# `BUILD_SHARED_LIBS=ON` build records the right soname without further edits.
set_target_properties(multifloats PROPERTIES
    VERSION   1.0.0
    SOVERSION 1
    # Hide all C++ symbols; only extern "C" dd_* functions are exported.
    CXX_VISIBILITY_PRESET    hidden
    VISIBILITY_INLINES_HIDDEN ON)

# Strip leaked inline-copy C++ symbols from the static archive.
# The compiler may emit out-of-line copies of inline header functions
# with global linkage; we localize them so only dd_* remain exported.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/localize_symbols.sh.in
               ${CMAKE_CURRENT_BINARY_DIR}/localize_symbols.sh @ONLY)
add_custom_command(TARGET multifloats POST_BUILD
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/localize_symbols.sh $<TARGET_FILE:multifloats>
    COMMENT "Localizing C++ symbols in libmultifloats (keeping dd_* only)")
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multifloats PRIVATE -O3 -flto -ffat-lto-objects)
    else()
        target_compile_options(multifloats PRIVATE -O3 -flto)
    endif()
else()
    target_compile_options(multifloats PRIVATE -O2)
endif()

# Install: library + header + export targets + Config.cmake
install(TARGETS multifloats
    EXPORT multifloatsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES multifloats.hh multifloats_c.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT multifloatsTargets
    FILE multifloatsTargets.cmake
    NAMESPACE multifloats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)

include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
"include(\"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets.cmake\")
")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)
