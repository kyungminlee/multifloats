cmake_minimum_required(VERSION 3.27)
project(multifloat LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(src)
# add_executable(test test/test.cc)

add_library(multifloat-fortran ./fsrc/multifloat.f90)
set_target_properties(multifloat-fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(multifloat-fortran PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

#add_executable(test-fortran ./fsrc/test.f90)
#target_link_libraries(test-fortran PRIVATE multifloat-fortran)

add_library(blas-multifloat
    ./blas/ddgemm.f90
    ./blas/ddtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloat-fortran)

add_executable(multifloat_fuzz ./fsrc/multifloat_fuzz.f90)
target_link_libraries(multifloat_fuzz PRIVATE multifloat-fortran)
add_executable(multifloat_test ./fsrc/multifloat_test.f90)
target_link_libraries(multifloat_test PRIVATE multifloat-fortran)
