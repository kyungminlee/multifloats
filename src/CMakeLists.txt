project(multifloats LANGUAGES CXX)

include(GNUInstallDirs)

set(MULTIFLOATS_PUBLIC_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

add_library(multifloats STATIC multifloats_math.cc multifloats_io.cc)
target_compile_features(multifloats PUBLIC cxx_std_17)
target_include_directories(multifloats
    PUBLIC $<BUILD_INTERFACE:${MULTIFLOATS_PUBLIC_INCLUDE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# ABI bookkeeping — parse MULTIFLOATS_ABI_VERSION straight out of the header
# so a bump there is the single source of truth. Without this, SOVERSION and
# the header constant drift silently and a future shared-library build would
# ship a wrongly-tagged soname. SOVERSION is only honored on shared libraries,
# but we set it anyway so `BUILD_SHARED_LIBS=ON` records the right soname.
file(STRINGS "${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats.h"
     _mf_abi_line REGEX "^#define[ \t]+MULTIFLOATS_ABI_VERSION[ \t]+[0-9]+")
if(NOT _mf_abi_line)
    message(FATAL_ERROR "multifloats: could not find MULTIFLOATS_ABI_VERSION in include/multifloats.h")
endif()
string(REGEX REPLACE "^#define[ \t]+MULTIFLOATS_ABI_VERSION[ \t]+([0-9]+).*$"
       "\\1" MULTIFLOATS_ABI_VERSION "${_mf_abi_line}")
set_target_properties(multifloats PROPERTIES
    VERSION   ${MULTIFLOATS_ABI_VERSION}.0.0
    SOVERSION ${MULTIFLOATS_ABI_VERSION}
    # Hide all C++ symbols; only extern "C" dd_* functions are exported.
    CXX_VISIBILITY_PRESET    hidden
    VISIBILITY_INLINES_HIDDEN ON)

# Strip non-API symbols from the static archive. The compiler emits
# out-of-line copies of some inline header functions + TU-internal helpers
# with (private-)external linkage; without this step they'd be visible as
# `T` entries in `nm libmultifloats.a` and resolvable by accident. The
# script's keep-list comes from two sources: the `MULTIFLOATS_API` scan of
# multifloats.h (C ABI) and a demangling scan of the archive itself for
# public `multifloats::NAME(...)` and `std::NAME<multifloats::float64x2>(...)`
# symbols. Opt out via `-DMULTIFLOATS_STRIP_SYMBOLS=OFF` for debugging.
option(MULTIFLOATS_STRIP_SYMBOLS
       "Localize non-API C++ symbols from libmultifloats.a" ON)
if(MULTIFLOATS_STRIP_SYMBOLS)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/localize_symbols.sh.in
                   ${CMAKE_CURRENT_BINARY_DIR}/localize_symbols.sh @ONLY)
    add_custom_command(TARGET multifloats POST_BUILD
        COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/localize_symbols.sh $<TARGET_FILE:multifloats>
        COMMENT "Localizing non-API symbols in libmultifloats")
endif()
# -O3 only — no -flto. LTO embeds compiler-specific bytecode (GIMPLE
# for GCC, opaque-format BC for Clang) into the static archive, which
# forces every downstream consumer onto a bytecode-compatible toolchain
# and breaks cross-compiler / cross-version linking against the
# extern "C" `*dd` ABI surface. The Fortran library keeps -flto
# (fsrc/CMakeLists.txt) because gfortran consumers already need a
# matching .mod / hidden-pointer ABI from the same gfortran build, so
# LTO adds no new constraint there. See doc/BENCHMARK.md "build" notes.
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    target_compile_options(multifloats PRIVATE -O3)
else()
    target_compile_options(multifloats PRIVATE -O2)
endif()

# Install: library + header + export targets + Config.cmake
install(TARGETS multifloats
    EXPORT multifloatsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES
    ${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats.h
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
