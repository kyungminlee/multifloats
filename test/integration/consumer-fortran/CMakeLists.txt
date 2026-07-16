# Minimal consumer project that exercises the installed multifloatsf
# package via find_package(). Drive from scripts/check_consumer_fortran.sh,
# which:
#   1. cmake --install the in-tree build into a scratch prefix.
#   2. Configure this project pointing CMAKE_PREFIX_PATH at that prefix.
#   3. Build, run, diff output against expected.
#
# The point is to prove that a fresh CMake project, with no knowledge of
# the multifloats build tree, can find the Fortran module and use both
# the bind(c) backend (through an intrinsic call) and an inline elemental
# operator. Output lines are compared verbatim by the driver script.

cmake_minimum_required(VERSION 3.18)
# Enable C and CXX alongside Fortran: the Fortran module's bind(c)
# interfaces link against libmultifloats (a C++-compiled library with a
# tagged-archive name under LTO), and multifloatsConfig.cmake selects the
# right archive by probing CMAKE_CXX_COMPILER_ID / _VERSION. With only
# Fortran enabled, the C library's selector sees an empty CXX tag and
# fails to match any installed build.
project(multifloats_fortran_consumer LANGUAGES CXX Fortran)

find_package(multifloatsf REQUIRED)

add_executable(consumer main.f90)
target_link_libraries(consumer PRIVATE multifloats::fortran_module)
