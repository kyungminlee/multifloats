project(multifloats LANGUAGES CXX)

include(GNUInstallDirs)

set(MULTIFLOATS_PUBLIC_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../include)

# ABI bookkeeping — parse MULTIFLOATS_ABI_VERSION straight out of the header
# so a bump there is the single source of truth. Without this, SOVERSION and
# the header constant drift silently and the shared-library build would ship a
# wrongly-tagged soname. SOVERSION is only honored on shared libraries, but the
# static archive carries it too so `BUILD_SHARED_LIBS=ON` records the right one.
file(STRINGS "${MULTIFLOATS_PUBLIC_INCLUDE_DIR}/multifloats/float64x2.h"
     _mf_abi_line REGEX "^#define[ \t]+MULTIFLOATS_ABI_VERSION[ \t]+[0-9]+")
if(NOT _mf_abi_line)
    message(FATAL_ERROR "multifloats: could not find MULTIFLOATS_ABI_VERSION in include/multifloats/float64x2.h")
endif()
string(REGEX REPLACE "^#define[ \t]+MULTIFLOATS_ABI_VERSION[ \t]+([0-9]+).*$"
       "\\1" MULTIFLOATS_ABI_VERSION "${_mf_abi_line}")

# Hide all C++ symbols; only extern "C" dd_* functions are exported.
# Optional because some consumers want full visibility for debugging,
# profiling tools that resolve inline-helper addresses, or when embedding
# the archive into a library that re-exports a wider surface.
option(MULTIFLOATS_HIDDEN_VISIBILITY
       "Hide non-API C++ symbols (-fvisibility=hidden) in libmultifloats" ON)

# Runtime FMA dispatch for the matmul kernels (see multifloats_math.cc). ON: the
# matmul implementation is compiled twice — a portable baseline and an AVX2+FMA
# copy whose panels vectorize to packed FMA (vfmadd…pd) — and the public entries
# pick one by CPUID at runtime, so a single binary uses packed FMA where present
# and runs on any CPU otherwise. OFF: a single baseline build (no second copy,
# no runtime branch). Effective only on x86-64 GCC — only GCC's vectorizer packs
# this compensated DD loop; a no-op on Clang/icx (they use the baseline), on
# AArch64 (FMA is baseline), and on Apple/Mach-O.
option(MULTIFLOATS_MM_DISPATCH
       "Runtime FMA-dispatch the matmul kernels (x86-64 GCC)" ON)

# Also build a versioned shared library (libmultifloats.so.<ABI> / .dylib)
# alongside the static archive, from the same compiled objects. OFF by default
# so an ordinary build produces only the static archive; the release turns it
# on in its non-LTO jobs to ship one portable shared object per platform.
option(MULTIFLOATS_BUILD_SHARED
       "Also build a versioned shared library (libmultifloats.so.N / .dylib)" OFF)

# NOTE on hardware FMA: the library deliberately does NOT compile with -mfma
# (or any FMA-enabling -march). Every DD multiply bottoms out in std::fma,
# which the compiler lowers to a call to glibc's `fma` — and that symbol is a
# GNU IFUNC: its resolver runs CPUID at load time and binds to the hardware
# FMA (vfmadd) implementation on capable CPUs, or a software fallback
# otherwise. So a single portable binary already uses hardware FMA at runtime
# where present and NEVER executes an illegal instruction on an older CPU.
# Statically inlining vfmadd (via -mfma) would only save the call/ifunc
# indirection (~12% on the transcendental kernels) but would SIGILL on any CPU
# without FMA — unacceptable for a library that must run anywhere. A consumer
# who fully controls their target CPUs can still pass -mfma/-march via
# CMAKE_CXX_FLAGS at their own risk.

# _mf_configure_objects(<target>): apply every setting that shapes the compiled
# objects. Called for BOTH the static archive and the optional shared library
# so their code is byte-identical (same fp semantics, visibility, and matmul
# dispatch) — this is the single source of truth for object-level config, so a
# new option added here automatically reaches both. Purely link-level and
# packaging differences (soname, hardening flags, LTO archiver) stay on the
# individual targets below.
function(_mf_configure_objects tgt)
    target_compile_features(${tgt} PUBLIC cxx_std_17)
    target_include_directories(${tgt}
        PUBLIC $<BUILD_INTERFACE:${MULTIFLOATS_PUBLIC_INCLUDE_DIR}>
               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    # PIC unconditionally so consumers can embed the static archive into
    # a downstream shared library (the `.so` linker requires every input
    # .o to be PIC) and so our own shared library can link. Negligible cost
    # on x86_64 / arm64.
    set_target_properties(${tgt} PROPERTIES POSITION_INDEPENDENT_CODE ON)

    if(MULTIFLOATS_HIDDEN_VISIBILITY)
        set_target_properties(${tgt} PROPERTIES
            CXX_VISIBILITY_PRESET     hidden
            VISIBILITY_INLINES_HIDDEN ON)
    endif()

    # Correctness — IntelLLVM. icpx defaults to -fp-model=fast at -O2+, which
    # reassociates the DD error-free transformations (2Sum / Fast2Sum / Dekker)
    # and contracts add-of-mul into FMAs — either corrupts the `lo`
    # compensation and breaks every DD operation. Unconditional (not gated on
    # build type): the optimizer rearranges expressions at -O1 too.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(${tgt} PRIVATE -fp-model=precise)
    endif()

    # Correctness — GCC/Clang. The default -ffp-contract=fast fuses an
    # incidental `a*b + c` into a single-rounded FMA on any FMA-capable target.
    # That happens with NO special flag on AArch64 (FMA is in the base ISA) and
    # whenever a consumer configures with -march=…/-mfma — silently changing the
    # rounding the error-free transformations are built around (a future
    # load-bearing `a*b+c` would be corrupted outright). Pin contraction off; it
    # comes after any inherited CMAKE_CXX_FLAGS on the command line, so it
    # neutralizes even -march=native. The explicit std::fma in two_prod is
    # unaffected (it is an fma by name, still lowered to glibc's ifunc-dispatched
    # hardware FMA), so no speed is lost — measured equal-or-faster. This is the
    # GCC/Clang analog of the IntelLLVM -fp-model=precise above, and makes every
    # build bit-identical to the reference regardless of target ISA.
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${tgt} PRIVATE -ffp-contract=off)
    endif()

    if(NOT MULTIFLOATS_MM_DISPATCH)
        target_compile_definitions(${tgt} PRIVATE MULTIFLOATS_NO_MM_DISPATCH)
    endif()
endfunction()

# ===========================================================================
# Static archive — the default artifact.
# ===========================================================================
add_library(multifloats STATIC multifloats_math.cc multifloats_io.cc)
_mf_configure_objects(multifloats)
set_target_properties(multifloats PROPERTIES
    VERSION   ${MULTIFLOATS_ABI_VERSION}.0.0
    SOVERSION ${MULTIFLOATS_ABI_VERSION})

# Namespace hygiene is enforced NON-destructively: the `check_exported_symbols`
# ctest (see the top-level CMakeLists.txt) asserts that the archive's only
# strong-global (T) symbols are the declared C ABI, the multifloats:: /
# multifloats::detail:: C++ API, and std::NAME<float64x2> specializations —
# failing the build (naming the offender) if a TU-internal helper escapes with
# external linkage. This replaced an earlier post-build strip (objcopy
# --keep-global-symbols + strip --strip-unneeded) that mutated the archive:
# that strip removed zero real leaks yet localized the weak COMDAT namespace
# copies, which GNU ld then discarded ("relocation refers to discarded
# section") — forcing STRIP_SYMBOLS=OFF on every non-LTO release. The audit
# keeps the guarantee without touching the archive, so it works identically
# for LTO and non-LTO builds. The primary source-level defense is unchanged:
# TU-internal helpers live in anonymous namespaces / `static` / detail::.

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
#
# Apple Clang is excluded explicitly: its `-ffat-lto-objects` probe
# returns success (the driver accepts the flag silently in trivial
# try_compile contexts) but the real compile errors with
# "unsupported option '-ffat-lto-objects' for target 'arm64-apple-darwin*'".
# The probe-vs-actual mismatch isn't worth working around — Apple Clang
# users get a non-fat LTO archive (consumable by Apple Clang at the same
# major version).
if(MULTIFLOATS_LTO_ON AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
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
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    string(REGEX MATCH "^([0-9]+)" _mf_cxx_ver "${CMAKE_CXX_COMPILER_VERSION}")
    set(_mf_cxx_tag "clang-${_mf_cxx_ver}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    # Apple Clang's version numbers are independent of upstream LLVM
    # Clang's, so tag them distinctly. (Apple Clang 16 ≠ LLVM Clang 16.)
    string(REGEX MATCH "^([0-9]+)" _mf_cxx_ver "${CMAKE_CXX_COMPILER_VERSION}")
    set(_mf_cxx_tag "appleclang-${_mf_cxx_ver}")
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
        OUTPUT_NAME "multifloats-lto-${MULTIFLOATS_CXX_TAG}")
    message(STATUS
        "multifloats: LTO on — archive will be libmultifloats-lto-${MULTIFLOATS_CXX_TAG}.a")
else()
    # Non-LTO archive is plain machine code, portable across compilers of the
    # same platform ABI, so it carries no compiler tag — only the `-nolto-`
    # axis label that keeps it distinct from the per-compiler LTO archives
    # when both are installed into one prefix.
    set_target_properties(multifloats PROPERTIES
        OUTPUT_NAME "multifloats-nolto")
    message(STATUS
        "multifloats: LTO off — archive will be libmultifloats-nolto.a")
endif()

# ===========================================================================
# Optional shared library — same compiled objects as the static archive,
# built as a versioned .so/.dylib.
#
# Its dynamic symbol table exports ONLY the extern "C" C ABI: hidden
# visibility keeps the multifloats:: C++ API and the weak COMDAT header
# inlines out of the dynsym, so the DSO surface is exactly the stable C ABI
# (SOVERSION == MULTIFLOATS_ABI_VERSION). One shared object per platform is
# built by the release's non-LTO jobs (plain machine code, compiler-agnostic).
# The `check_shared_exports` ctest audits the resulting dynsym.
# ===========================================================================
if(MULTIFLOATS_BUILD_SHARED)
    add_library(multifloats_shared SHARED multifloats_math.cc multifloats_io.cc)
    _mf_configure_objects(multifloats_shared)
    set_target_properties(multifloats_shared PROPERTIES
        OUTPUT_NAME multifloats
        VERSION   ${MULTIFLOATS_ABI_VERSION}.0.0
        SOVERSION ${MULTIFLOATS_ABI_VERSION})

    # Match the static archive's LTO setting when both are requested together.
    # (In the release the shared library is built only in the non-LTO jobs, so
    # this is normally a no-op; it keeps a hand `-DMULTIFLOATS_BUILD_SHARED=ON
    # -DMULTIFLOATS_USE_LTO=ON` build consistent.)
    if(MULTIFLOATS_LTO_ON)
        set_target_properties(multifloats_shared PROPERTIES
            INTERPROCEDURAL_OPTIMIZATION_RELEASE        TRUE
            INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE
            INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL     TRUE)
    endif()

    # Link-time hardening + a fully self-contained shared object:
    #   ELF   : -z relro,-z now → full RELRO + immediate binding (BIND_NOW), so
    #           the GOT is read-only after load and there is no lazy-PLT
    #           surface; --no-undefined fails the link if any symbol is
    #           unresolved, forcing libstdc++/libm to be recorded NEEDED.
    #   Mach-O: -bind_at_load is the -z now analog (bind everything at load);
    #           ld64 already errors on undefined symbols by default.
    if(APPLE)
        target_link_options(multifloats_shared PRIVATE LINKER:-bind_at_load)
    elseif(UNIX)
        target_link_options(multifloats_shared PRIVATE
            LINKER:-z,relro LINKER:-z,now LINKER:--no-undefined)
        target_link_libraries(multifloats_shared PRIVATE m)
    endif()

    # Installed into libdir with the usual versioned soname + dev symlinks
    # (libmultifloats.so.N.0.0 ← .so.N ← .so). Deliberately NOT added to the
    # multifloatsTargets export set: find_package(multifloats) keeps routing
    # to the static archive via its LTO/nolto selector, and the shared object
    # is for direct `-lmultifloats` / dlopen binary consumers.
    install(TARGETS multifloats_shared
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    message(STATUS
        "multifloats: shared library on — will build libmultifloats.so.${MULTIFLOATS_ABI_VERSION} / .dylib")
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

# Targets file: LTO builds are per-compiler-tagged (`-lto-<tag>`) so several
# can coexist under one prefix; the non-LTO build carries the compiler-agnostic
# `-nolto` label (its archive is plain machine code, portable across compilers).
if(MULTIFLOATS_LTO_ON)
    set(_mf_targets_filename "multifloatsTargets-lto-${MULTIFLOATS_CXX_TAG}.cmake")
else()
    set(_mf_targets_filename "multifloatsTargets-nolto.cmake")
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
#   2. Looks for an exact-match `multifloatsTargets-lto-<tag>.cmake` (an LTO
#      build keyed to that exact compiler).
#   3. Falls back to `multifloatsTargets-nolto.cmake` (the compiler-agnostic
#      non-LTO build, which contains plain machine code and is portable).
#   4. Errors out only if neither exists.
include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/multifloatsConfig.cmake"
"# multifloatsConfig.cmake — auto-generated, best-fit selector.
# Tries an exact compiler-tag match first (LTO build), falls back to
# the compiler-agnostic non-LTO build.
${_MF_CXX_DETECT_CODE}
set(_mf_targets \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets-lto-\${_mf_cxx_tag}.cmake\")
if(NOT EXISTS \"\${_mf_targets}\")
    set(_mf_targets \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets-nolto.cmake\")
endif()
if(NOT EXISTS \"\${_mf_targets}\")
    file(GLOB _mf_available \"\${CMAKE_CURRENT_LIST_DIR}/multifloatsTargets-lto-*.cmake\")
    set(_mf_avail_tags \"\")
    foreach(_mf_f IN LISTS _mf_available)
        get_filename_component(_mf_fn \"\${_mf_f}\" NAME)
        # Strip per-config helpers (CMake generates <export>-<config>.cmake
        # alongside the main file; only the main file is a compiler tag).
        if(_mf_fn MATCHES \"-(noconfig|release|debug|relwithdebinfo|minsizerel)\\\\.cmake$\")
            continue()
        endif()
        string(REGEX REPLACE \"^multifloatsTargets-lto-(.*)\\\\.cmake$\" \"\\\\1\" _mf_t \"\${_mf_fn}\")
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
