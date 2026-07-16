# Standalone example project: two tiny programs that compute π to full
# double-double precision, one through the C++ API and one through the
# Fortran module. This is NOT part of the library build — configure it
# separately against an installed multifloats / multifloatsf package:
#
#   cmake --build build                                   # build the library
#   cmake --install build --prefix /path/to/prefix        # stage it
#   cmake -S example -B build/example -DCMAKE_PREFIX_PATH=/path/to/prefix
#   cmake --build build/example
#   ./build/example/pi_cxx ; ./build/example/pi_f
#
# CXX is enabled even for the Fortran example because
# multifloats::fortran_module links libmultifloats (a C++-compiled archive)
# and its Config selector probes the CXX compiler tag.
cmake_minimum_required(VERSION 3.18)
project(multifloats_examples LANGUAGES CXX Fortran)

find_package(multifloats  REQUIRED)   # C / C++ kernels
find_package(multifloatsf REQUIRED)   # Fortran module

add_executable(pi_cxx pi.cpp)
target_link_libraries(pi_cxx PRIVATE multifloats::multifloats)

add_executable(pi_f pi.f90)
target_link_libraries(pi_f PRIVATE multifloats::fortran_module)
