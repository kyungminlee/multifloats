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
file(STRINGS "${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats/float64x2.h"
     _mf_abi_line REGEX "^#define[ \t]+MULTIFLOATS_ABI_VERSION[ \t]+[0-9]+")
if(NOT _mf_abi_line)
    message(FATAL_ERROR "multifloats: could not find MULTIFLOATS_ABI_VERSION in include/multifloats/float64x2.h")
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
# Correctness — IntelLLVM. icpx defaults to -fp-model=fast at -O2+, which
# reassociates the DD error-free transformations (2Sum / Fast2Sum / Dekker)
# and contracts add-of-mul into FMAs — either corrupts the `lo`
# compensation and breaks every DD operation. Unconditional (not gated on
# build type): the optimizer rearranges expressions at -O1 too.
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(multifloats PRIVATE -fp-model=precise)
endif()

# LTO via CMake's standard abstraction, gated on the project-wide
# MULTIFLOATS_USE_LTO option (defined in the top-level CMakeLists.txt).
# Per-config so a Debug build leaves IPO off (LTO destroys debuggability
# and noticeably slows the link). check_ipo_supported errors out by
# default on toolchains that don't do IPO; LANGUAGES CXX scopes the
# probe to the language we care about.
set(MULTIFLOATS_LTO_ON FALSE)
if(MULTIFLOATS_USE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT MULTIFLOATS_IPO_SUPPORTED LANGUAGES CXX)
    if(MULTIFLOATS_IPO_SUPPORTED)
        set_target_properties(multifloats PROPERTIES
            INTERPROCEDURAL_OPTIMIZATION_RELEASE        TRUE
            INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE
            INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL     TRUE)
        set(MULTIFLOATS_LTO_ON TRUE)
    endif()
endif()

# Fat LTO: emit both LTO bitcode and ordinary object code so the archive
# remains linkable by a non-LTO linker or a different-vendor toolchain
# (e.g. a GCC-built libmultifloats.a consumed by a Clang-driven link).
# CMake doesn't model fat-LTO directly, so we probe + add the flag manually,
# scoped to the optimized configs that actually enable IPO.
if(MULTIFLOATS_LTO_ON)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-ffat-lto-objects MULTIFLOATS_CXX_HAS_FAT_LTO)
    if(MULTIFLOATS_CXX_HAS_FAT_LTO)
        target_compile_options(multifloats PRIVATE
            $<$<CONFIG:Release,RelWithDebInfo,MinSizeRel>:-ffat-lto-objects>)
    endif()
endif()

# Install: library + header + export targets + Config.cmake
install(TARGETS multifloats
    EXPORT multifloatsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES
    ${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES
    ${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats/float64x2.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/multifloats)
install(EXPORT multifloatsTargets
    FILE multifloatsTargets.cmake
    NAMESPACE multifloats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)

include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
"include(\"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets.cmake\")
")
# Version file so `find_package(multifloats <N> REQUIRED)` actually
# enforces the requested major. SameMajorVersion: any installed N.x.x
# satisfies a request for N (matching what MULTIFLOATS_ABI_VERSION ==
# SOVERSION encodes — bump on any C-ABI break).
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfigVersion.cmake"
    VERSION ${MULTIFLOATS_ABI_VERSION}.0.0
    COMPATIBILITY SameMajorVersion)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)
