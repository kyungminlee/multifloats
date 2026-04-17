cmake_minimum_required(VERSION 3.27)

# Pin the compilers to Homebrew GCC 15 by default. Apple Clang lacks
# __float128 / libquadmath, and Apple-shipped LLVM Flang on Apple Silicon
# does not enable REAL(KIND=16). Override on the command line with
# -DCMAKE_CXX_COMPILER=... / -DCMAKE_Fortran_COMPILER=...
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

project(multifloats-dev LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_program(FYPP_EXECUTABLE fypp REQUIRED)

# =============================================================================
# Library subprojects
# =============================================================================
add_subdirectory(src)   # project(multifloats)  → libmultifloats.a
add_subdirectory(fsrc)  # project(multifloatsf) → libmultifloatsf-<compiler>.a

# =============================================================================
# Tests and benchmarks (not installed)
# =============================================================================

# C++ tests — use __float128 / libquadmath as reference (GCC only).
add_executable(cpp_test ./test/test.cc)
target_compile_options(cpp_test PRIVATE -O3 -fext-numeric-literals)
target_link_libraries(cpp_test PRIVATE multifloats quadmath)
add_test(NAME cpp_test COMMAND cpp_test)

add_executable(cpp_fuzz ./test/fuzz.cc)
target_compile_options(cpp_fuzz PRIVATE -O3 -fext-numeric-literals)
target_link_libraries(cpp_fuzz PRIVATE multifloats quadmath)
add_test(NAME cpp_fuzz COMMAND cpp_fuzz)

add_executable(cpp_bench ./test/bench.cc)
target_compile_options(cpp_bench PRIVATE -O3 -fext-numeric-literals)
target_link_libraries(cpp_bench PRIVATE multifloats quadmath)

# Fortran tests
add_executable(fortran_precision ./test/precision.f90)
target_compile_options(fortran_precision PRIVATE -O3)
target_link_libraries(fortran_precision PRIVATE multifloatsf)
add_test(NAME fortran_precision COMMAND fortran_precision)

add_executable(fortran_test ./test/test.f90)
target_compile_options(fortran_test PRIVATE -O3)
target_link_libraries(fortran_test PRIVATE multifloatsf)
add_test(NAME fortran_test COMMAND fortran_test)

add_executable(fortran_fuzz ./test/fuzz.f90)
target_compile_options(fortran_fuzz PRIVATE -O3)
target_link_libraries(fortran_fuzz PRIVATE multifloatsf)
add_test(NAME fortran_fuzz COMMAND fortran_fuzz)

add_executable(fortran_bench ./test/bench.f90)
target_compile_options(fortran_bench PRIVATE -O3)
target_link_libraries(fortran_bench PRIVATE multifloatsf)

# BLAS extensions (not installed)
add_library(blas-multifloat STATIC
    ./blas/wgemm.f90
    ./blas/wtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloatsf)

# Fortran library uses C++ DD kernels via bind(C) (extern "C" in multifloats_math.cc)
target_link_libraries(multifloatsf PUBLIC multifloats)

add_executable(fortran_bench_c ./test/bench_c.f90)
target_compile_options(fortran_bench_c PRIVATE -O3)
target_link_libraries(fortran_bench_c PRIVATE multifloatsf)

add_library(dd-bindc OBJECT ./test/dd_bindc.f90)
target_compile_options(dd-bindc PRIVATE -O3)
set_target_properties(dd-bindc PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fmod/test)
target_include_directories(dd-bindc PUBLIC ${CMAKE_BINARY_DIR}/fmod/test)

add_executable(fortran_bench_abi ./test/bench_abi.f90)
target_compile_options(fortran_bench_abi PRIVATE -O3)
target_link_libraries(fortran_bench_abi PRIVATE multifloatsf dd-bindc)
