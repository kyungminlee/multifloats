cmake_minimum_required(VERSION 3.27)

# Pin the compilers to Homebrew GCC 15 by default. Apple Clang lacks
# __float128 / libquadmath, and Apple-shipped LLVM Flang on Apple Silicon
# does not enable REAL(KIND=16). Both are required by the multifloats
# Fortran module and the C++ test harness. Override on the command line
# with -DCMAKE_CXX_COMPILER=... / -DCMAKE_Fortran_COMPILER=...
if(NOT DEFINED CMAKE_CXX_COMPILER)
    find_program(_MF_GXX NAMES g++-15 g++-14 g++-13
        PATHS /opt/homebrew/bin /usr/local/bin /usr/bin)
    if(_MF_GXX)
        set(CMAKE_CXX_COMPILER ${_MF_GXX})
    endif()
endif()
if(NOT DEFINED CMAKE_Fortran_COMPILER)
    find_program(_MF_GFORTRAN NAMES gfortran-15 gfortran-14 gfortran-13
        PATHS /opt/homebrew/bin /usr/local/bin /usr/bin)
    if(_MF_GFORTRAN)
        set(CMAKE_Fortran_COMPILER ${_MF_GFORTRAN})
    endif()
endif()

project(multifloat LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

include_directories(src)
find_program(FYPP_EXECUTABLE fypp REQUIRED)

set(MULTIFLOATS_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(MULTIFLOATS_GENERATED_SRC ${MULTIFLOATS_GENERATED_DIR}/multifloats.f90)

# Detect REAL(16) support: gfortran has it, flang on ARM64 does not.
include(CheckFortranSourceCompiles)
check_fortran_source_compiles(
    "program t; real(16) :: x; x = 1.0; end program"
    FORTRAN_HAS_REAL16
    SRC_EXT f90)
set(_FYPP_DEFS "")
if(NOT FORTRAN_HAS_REAL16)
    list(APPEND _FYPP_DEFS "-DMULTIFLOATS_NO_QUAD")
endif()

add_custom_command(
    OUTPUT ${MULTIFLOATS_GENERATED_SRC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MULTIFLOATS_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE} ${_FYPP_DEFS} ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp ${MULTIFLOATS_GENERATED_SRC}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp
    VERBATIM)

add_library(multifloatsf STATIC ${MULTIFLOATS_GENERATED_SRC})
set_target_properties(multifloatsf PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(multifloatsf PUBLIC
    ${CMAKE_BINARY_DIR}/modules
    ${CMAKE_CURRENT_SOURCE_DIR}/src)
# Fat LTO for GCC: embed both GIMPLE IR and native code so ar preserves
# LTO in .a archives. Flang/LLVM uses -flto=full (thin LTO also works).
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multifloatsf PUBLIC -flto -ffat-lto-objects)
        target_link_options(multifloatsf PUBLIC -flto)
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
        target_compile_options(multifloatsf PUBLIC -flto)
        target_link_options(multifloatsf PUBLIC -flto)
    endif()
endif()

add_executable(fortran_precision ./test/precision.f90)
target_link_libraries(fortran_precision PRIVATE multifloatsf)
add_test(NAME fortran_precision COMMAND fortran_precision)

add_library(blas-multifloat
    ./blas/wgemm.f90
    ./blas/wtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloatsf)

add_executable(fortran_fuzz ./test/fuzz.f90)
target_link_libraries(fortran_fuzz PRIVATE multifloatsf)
add_test(NAME fortran_fuzz COMMAND fortran_fuzz)

add_executable(fortran_test ./test/test.f90)
target_link_libraries(fortran_test PRIVATE multifloatsf)
add_test(NAME fortran_test COMMAND fortran_test)

add_executable(fortran_bench ./test/bench.f90)
target_link_libraries(fortran_bench PRIVATE multifloatsf)

# DD math kernel implementations (separated from the header-only template
# library so that the large function bodies are compiled once).
add_library(multifloats STATIC ./src/multifloats_math.cc)
target_include_directories(multifloats PRIVATE src)
target_compile_options(multifloats PRIVATE -O3)
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multifloats PRIVATE -flto -ffat-lto-objects)
    else()
        target_compile_options(multifloats PRIVATE -flto)
    endif()
endif()

# C++ tests for src/multifloats.hh — use __float128 / libquadmath as a
# high-precision reference (provided by GCC, not Apple Clang).
add_executable(cpp_test ./test/test.cc)
target_compile_options(cpp_test PRIVATE -fext-numeric-literals)
target_link_libraries(cpp_test PRIVATE multifloats quadmath)
add_test(NAME cpp_test COMMAND cpp_test)

add_executable(cpp_fuzz ./test/fuzz.cc)
target_compile_options(cpp_fuzz PRIVATE -fext-numeric-literals)
target_link_libraries(cpp_fuzz PRIVATE multifloats quadmath)
add_test(NAME cpp_fuzz COMMAND cpp_fuzz)

add_executable(cpp_bench ./test/bench.cc)
target_compile_options(cpp_bench PRIVATE -fext-numeric-literals)
target_link_libraries(cpp_bench PRIVATE multifloats quadmath)

# C-ABI wrapper library (compiled as C++ with full inlining, exposed as
# extern "C"). Used by bench_c to compare gfortran's derived-type ABI
# against the platform C calling convention for the same DD kernels.
add_library(multifloat-c STATIC ./src/multifloats_c.cc)
target_include_directories(multifloat-c PRIVATE src)
target_compile_options(multifloat-c PRIVATE -O3)
target_link_libraries(multifloat-c PRIVATE multifloats)

add_executable(fortran_bench_c ./test/bench_c.f90)
target_link_libraries(fortran_bench_c PRIVATE multifloatsf multifloat-c)

add_library(dd-bindc OBJECT ./test/dd_bindc.f90)
target_compile_options(dd-bindc PRIVATE -O3)
set_target_properties(dd-bindc PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/test)

add_executable(fortran_bench_abi ./test/bench_abi.f90)
target_link_libraries(fortran_bench_abi PRIVATE multifloatsf multifloat-c dd-bindc)

# =============================================================================
# Install rules
# =============================================================================
include(GNUInstallDirs)

# C++ library + header
install(TARGETS multifloats
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES src/multifloats.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Fortran library + .mod file (versioned by compiler ID + major version,
# since .mod format is compiler-specific and version-specific).
install(TARGETS multifloatsf
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" _fc_id)
set(_fc_mod_dir "${CMAKE_INSTALL_LIBDIR}/multifloats/${_fc_id}-${CMAKE_Fortran_COMPILER_VERSION}")
install(FILES ${CMAKE_BINARY_DIR}/modules/multifloats.mod
    DESTINATION ${_fc_mod_dir})
