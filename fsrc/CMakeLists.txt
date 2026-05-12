project(multifloatsf LANGUAGES Fortran)

include(GNUInstallDirs)
include(FortranCompiler)

# Detect usable REAL(16) support on the *build* host. Governs which
# pre-expanded .f90 variant we use for the in-tree multifloatsf target
# (tests, bench). Consumers making a source-build via find_package()
# re-run the same probe themselves in multifloatsfConfig.cmake and pick
# whichever of the two shipped variants matches their compiler.
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("
    program t
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      real(16) :: x, y
      real(8)  :: d
      x = 1.0
      y = x + x
      d = real(x, 8)
      if (ieee_is_finite(x)) y = y * x
    end program
" FORTRAN_HAS_REAL16 SRC_EXT f90)

set(_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)

# Generate both variants unconditionally. Both go into the install image
# so cross-compilation and consumer-side compiler switches are safe; the
# build host's REAL(16) status only decides which variant the in-tree
# multifloatsf target compiles.
set(_SRC_QUAD   ${_GENERATED_DIR}/multifloats-quad.f90)
set(_SRC_NOQUAD ${_GENERATED_DIR}/multifloats-noquad.f90)

add_custom_command(
    OUTPUT ${_SRC_QUAD}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
        ${_SRC_QUAD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
    VERBATIM)
add_custom_command(
    OUTPUT ${_SRC_NOQUAD}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE} -DMULTIFLOATS_NO_QUAD
        ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
        ${_SRC_NOQUAD}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
    VERBATIM)

# Aggregate target so both variants are built whenever anything in this
# tree depends on "the fypp output" (needed before install). The `ALL`
# attachment ensures both land during plain `cmake --build`, but CI
# release jobs that invoke `cmake --build build --target multifloats
# --target multifloatsf` (LTO non-GCC matrix entries) bypass ALL-linked
# targets — wire the aggregate to `multifloatsf` explicitly as well so
# both .f90 variants exist by the time `cmake --install` runs.
add_custom_target(multifloatsf_sources ALL
    DEPENDS ${_SRC_QUAD} ${_SRC_NOQUAD})

# In-tree multifloatsf library: picks the variant that matches the build
# host's REAL(16) support. scripts/check_fortran_abi_sync.sh and the
# cross-check tests point at this concrete file — keep a stable
# `generated/multifloats.f90` symlink-style dispatch for them.
if(FORTRAN_HAS_REAL16)
    set(_GENERATED_SRC ${_SRC_QUAD})
else()
    set(_GENERATED_SRC ${_SRC_NOQUAD})
endif()
# Compatibility alias expected by scripts/check_fortran_abi_sync.sh which
# greps generated/multifloats.f90 for the bind(c) surface.
set(_COMPAT_SRC ${_GENERATED_DIR}/multifloats.f90)
add_custom_command(
    OUTPUT ${_COMPAT_SRC}
    COMMAND ${CMAKE_COMMAND} -E copy ${_GENERATED_SRC} ${_COMPAT_SRC}
    DEPENDS ${_GENERATED_SRC}
    VERBATIM)
add_custom_target(multifloatsf_compat_alias ALL DEPENDS ${_COMPAT_SRC})

add_library(multifloatsf STATIC ${_GENERATED_SRC})
add_dependencies(multifloatsf multifloatsf_sources)
fortran_module_layout(multifloatsf)
target_include_directories(multifloatsf PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# PIC unconditionally so consumers can embed the static archive into a
# downstream shared library (matches src/CMakeLists.txt; see rationale
# there). Negligible cost on x86_64 / arm64.
set_target_properties(multifloatsf PROPERTIES
    POSITION_INDEPENDENT_CODE ON)

# Correctness — IntelLLVM. ifx defaults to -fp-model=fast at -O2+. DD
# error-free transformations require strict IEEE associativity and no
# FMA contraction — see src/CMakeLists.txt for the full rationale.
# Unconditional (not gated on build type): the optimizer rearranges at
# -O1 too.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(multifloatsf PRIVATE -fp-model=precise)
endif()

# LTO via CMake's IPO abstraction, gated on the project-wide
# MULTIFLOATS_USE_LTO option (defined in the top-level CMakeLists.txt)
# AND on GNU / LLVMFlang because other Fortran compilers (Intel ifx,
# NVIDIA, NAG, Cray) either don't have a usable LTO mode or produce LTO
# bytecode no other compiler can consume, defeating the cross-toolchain
# ABI surface this library exposes. The CheckIPOSupported probe alone
# wouldn't catch the interop concern.
set(MULTIFLOATSF_LTO_ON FALSE)
if(MULTIFLOATS_USE_LTO
   AND (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
        OR CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang"))
    include(CheckIPOSupported)
    check_ipo_supported(RESULT MULTIFLOATSF_IPO_SUPPORTED LANGUAGES Fortran)
    if(MULTIFLOATSF_IPO_SUPPORTED)
        set_target_properties(multifloatsf PROPERTIES
            INTERPROCEDURAL_OPTIMIZATION_RELEASE        TRUE
            INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE
            INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL     TRUE)
        set(MULTIFLOATSF_LTO_ON TRUE)
    endif()
endif()

# Fat LTO: keep regular object code alongside the LTO bitcode so the
# installed archive can be linked by a different compiler/toolchain than
# the one that built it (non-LTO or LTO-incompatible linker driver
# falls back to the embedded native objects). CMake doesn't model
# fat-LTO directly, so probe + add manually, scoped to the optimized
# configs that actually enable IPO.
if(MULTIFLOATSF_LTO_ON)
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag(-ffat-lto-objects MULTIFLOATSF_HAS_FAT_LTO)
    if(MULTIFLOATSF_HAS_FAT_LTO)
        target_compile_options(multifloatsf PRIVATE
            $<$<CONFIG:Release,RelWithDebInfo,MinSizeRel>:-ffat-lto-objects>)
    endif()
endif()

# =============================================================================
# Install
# =============================================================================
# Default distribution model is source: the two fypp-expanded .f90 variants
# plus a Config.cmake that picks and compiles at consumer find_package()
# time. Precompiled-mod install is opt-in (see top-level
# MULTIFLOATSF_INSTALL_PRECOMPILED_MOD option).

# Where the .f90 lives in the install tree. Under share/ rather than lib/
# because it's architecture-independent and not a linkable artifact.
set(MULTIFLOATSF_SHARE_DIR share/multifloatsf/src
    CACHE STRING "Install subdir (below <prefix>) for the fypp-expanded Fortran sources")

install(FILES ${_SRC_QUAD} ${_SRC_NOQUAD}
    DESTINATION ${MULTIFLOATSF_SHARE_DIR})
# Also ship the fypp template itself for consumers who want to regenerate
# (e.g. custom fypp macros). Optional — most consumers use the
# pre-expanded .f90s and never touch fypp.
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
    DESTINATION ${MULTIFLOATSF_SHARE_DIR})

# Opt-in: also install a precompiled compiler-tagged archive + .mod for
# consumers whose Fortran compiler matches a shipped tag. Config.cmake
# prefers this path when available and falls back to source-build
# otherwise, so the two modes compose cleanly. We cannot call
# fortran_install_library() here because that helper also writes a
# multifloatsfConfig.cmake of its own (best-fit-or-fail, no source
# fallback); it would clobber the unified Config we install below. So
# we inline just the tagged-archive + tagged-export install pieces.
if(MULTIFLOATSF_INSTALL_PRECOMPILED_MOD)
    set_target_properties(multifloatsf PROPERTIES
        OUTPUT_NAME "multifloatsf-${FORTRAN_COMPILER_TAG}")
    install(TARGETS multifloatsf
        EXPORT multifloatsfTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(EXPORT multifloatsfTargets
        FILE multifloatsfTargets-${FORTRAN_COMPILER_TAG}.cmake
        NAMESPACE multifloats::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloatsf)
    fortran_install_modules(multifloatsf)
endif()

# Config.cmake: unified selector (precompiled-first, source-fallback).
# The `@FC_DETECT_CODE@` placeholder in the template is expanded from
# the _FC_DETECT_CODE snippet that FortranCompiler.cmake wrote at its
# own configure-time — so producer and consumer compute identical tags.
include(CMakePackageConfigHelpers)
set(FC_DETECT_CODE "${_FC_DETECT_CODE}")
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/../cmake/multifloatsfConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/multifloatsfConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloatsf
    PATH_VARS MULTIFLOATSF_SHARE_DIR)
unset(FC_DETECT_CODE)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/multifloatsfConfigVersion.cmake
    VERSION ${MULTIFLOATS_ABI_VERSION}.0.0
    COMPATIBILITY SameMajorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/multifloatsfConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/multifloatsfConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multifloatsf)
