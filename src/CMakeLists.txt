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
    VISIBILITY_INLINES_HIDDEN ON
    # PIC unconditionally so consumers can embed the static archive into
    # a downstream shared library (the `.so` linker requires every input
    # .o to be PIC). Negligible cost on x86_64 / arm64.
    POSITION_INDEPENDENT_CODE ON)

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

# When LTO is on, the static archive embeds compiler-version-locked
# bytecode (GCC GIMPLE / Clang BC). Tag the OUTPUT_NAME so two builds
# from different compilers — or different majors of the same compiler —
# can install side by side under the same prefix without overwriting
# each other. Mirrors the Fortran library's `libmultifloatsf-gfortran-15.a`
# pattern (see cmake/FortranCompiler.cmake), on a different ABI axis
# (LTO bytecode here, .mod format there).
#
# The detection snippet is shared between producer (this file, at build
# time) and consumer (the generated multifloatsConfig.cmake, at
# find_package time) so both compute identical tags. Without that, the
# consumer would derive a different name than the installed archive
# carries and find_package would silently miss the build.
set(_MF_CXX_DETECT_CODE [[
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    string(REGEX MATCH "^([0-9]+)" _mf_cxx_ver "${CMAKE_CXX_COMPILER_VERSION}")
    set(_mf_cxx_tag "gcc-${_mf_cxx_ver}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
       OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    string(REGEX MATCH "^([0-9]+)" _mf_cxx_ver "${CMAKE_CXX_COMPILER_VERSION}")
    set(_mf_cxx_tag "clang-${_mf_cxx_ver}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    string(REGEX MATCH "^([0-9]+\\.[0-9]+)" _mf_cxx_ver "${CMAKE_CXX_COMPILER_VERSION}")
    set(_mf_cxx_tag "intel-${_mf_cxx_ver}")
else()
    string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" _mf_cxx_id)
    set(_mf_cxx_tag "${_mf_cxx_id}-${CMAKE_CXX_COMPILER_VERSION}")
    unset(_mf_cxx_id)
endif()
unset(_mf_cxx_ver)
]])

if(MULTIFLOATS_LTO_ON)
    cmake_language(EVAL CODE "${_MF_CXX_DETECT_CODE}")
    set(MULTIFLOATS_CXX_TAG "${_mf_cxx_tag}")
    unset(_mf_cxx_tag)
    set_target_properties(multifloats PROPERTIES
        OUTPUT_NAME "multifloats-${MULTIFLOATS_CXX_TAG}")
    message(STATUS
        "multifloats: LTO on — archive will be libmultifloats-${MULTIFLOATS_CXX_TAG}.a")
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

# Targets file: tagged when LTO is on so multiple per-compiler builds
# can coexist under the same install prefix. Untagged otherwise (a
# non-LTO archive contains only plain machine code and is portable).
if(MULTIFLOATS_LTO_ON)
    set(_mf_targets_filename "multifloatsTargets-${MULTIFLOATS_CXX_TAG}.cmake")
else()
    set(_mf_targets_filename "multifloatsTargets.cmake")
endif()
install(EXPORT multifloatsTargets
    FILE ${_mf_targets_filename}
    NAMESPACE multifloats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloats)

# Config.cmake: best-fit selector. Always generated the same way (whether
# this build is LTO-on or LTO-off) so multiple builds installed into the
# same prefix don't overwrite each other's Config with a less-capable
# variant. The selector at consumer find_package() time:
#   1. Derives the consumer's compiler tag (gcc-15 / clang-19 / ...).
#   2. Looks for an exact-match `multifloatsTargets-<tag>.cmake` (an LTO
#      build keyed to that exact compiler).
#   3. Falls back to `multifloatsTargets.cmake` (an untagged non-LTO
#      build, which contains plain machine code and is portable).
#   4. Errors out only if neither exists.
include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
"# multifloatsConfig.cmake — auto-generated, best-fit selector.
# Tries an exact compiler-tag match first (LTO build), falls back to
# the untagged variant (portable non-LTO build).
${_MF_CXX_DETECT_CODE}
set(_mf_targets \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets-\${_mf_cxx_tag}.cmake\")
if(NOT EXISTS \"\${_mf_targets}\")
    set(_mf_targets \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets.cmake\")
endif()
if(NOT EXISTS \"\${_mf_targets}\")
    file(GLOB _mf_available \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets-*.cmake\")
    set(_mf_avail_tags \"\")
    foreach(_mf_f IN LISTS _mf_available)
        get_filename_component(_mf_fn \"\${_mf_f}\" NAME)
        # Strip per-config helpers (CMake generates <export>-<config>.cmake
        # alongside the main file; only the main file is a compiler tag).
        if(_mf_fn MATCHES \"-(noconfig|release|debug|relwithdebinfo|minsizerel)\\\\.cmake$\")
            continue()
        endif()
        string(REGEX REPLACE \"^multifloatsTargets-(.*)\\\\.cmake$\" \"\\\\1\" _mf_t \"\${_mf_fn}\")
        list(APPEND _mf_avail_tags \"\${_mf_t}\")
    endforeach()
    list(SORT _mf_avail_tags)
    list(JOIN _mf_avail_tags \", \" _mf_avail_list)
    set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
    set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
        \"multifloats: no build matching compiler '\${_mf_cxx_tag}' and no portable non-LTO build installed. LTO tags available: [\${_mf_avail_list}]. Rebuild with the same compiler, or install a non-LTO build (-DMULTIFLOATS_USE_LTO=OFF).\")
    unset(_mf_cxx_tag)
    unset(_mf_targets)
    unset(_mf_available)
    unset(_mf_avail_tags)
    unset(_mf_avail_list)
    return()
endif()
include(\"\${_mf_targets}\")
unset(_mf_cxx_tag)
unset(_mf_targets)
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
