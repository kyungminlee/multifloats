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

add_custom_command(
    OUTPUT ${MULTIFLOATS_GENERATED_SRC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MULTIFLOATS_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp ${MULTIFLOATS_GENERATED_SRC}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp
    VERBATIM)

add_library(multifloat-fortran ${MULTIFLOATS_GENERATED_SRC})
set_target_properties(multifloat-fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(multifloat-fortran PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

add_executable(fortran_precision ./test/precision.f90)
target_link_libraries(fortran_precision PRIVATE multifloat-fortran)
add_test(NAME fortran_precision COMMAND fortran_precision)

add_library(blas-multifloat
    ./blas/wgemm.f90
    ./blas/wtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloat-fortran)

add_executable(fortran_fuzz ./test/fuzz.f90)
target_link_libraries(fortran_fuzz PRIVATE multifloat-fortran)
add_test(NAME fortran_fuzz COMMAND fortran_fuzz)

add_executable(fortran_test ./test/test.f90)
target_link_libraries(fortran_test PRIVATE multifloat-fortran)
add_test(NAME fortran_test COMMAND fortran_test)

# C++ tests for src/multifloats.hh — use __float128 / libquadmath as a
# high-precision reference (provided by GCC, not Apple Clang).
add_executable(cpp_test ./test/test.cc)
target_compile_options(cpp_test PRIVATE -fext-numeric-literals)
target_link_libraries(cpp_test PRIVATE quadmath)
add_test(NAME cpp_test COMMAND cpp_test)

add_executable(cpp_fuzz ./test/fuzz.cc)
target_compile_options(cpp_fuzz PRIVATE -fext-numeric-literals)
target_link_libraries(cpp_fuzz PRIVATE quadmath)
add_test(NAME cpp_fuzz COMMAND cpp_fuzz)
