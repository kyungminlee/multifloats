cmake_minimum_required(VERSION 3.27)
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

add_executable(test-fortran ./fsrc/test.f90)
target_link_libraries(test-fortran PRIVATE multifloat-fortran)
add_test(NAME test-fortran COMMAND test-fortran)

add_library(blas-multifloat
    ./blas/ddgemm.f90
    ./blas/ddtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloat-fortran)

add_executable(multifloat_fuzz ./fsrc/multifloat_fuzz.f90)
target_link_libraries(multifloat_fuzz PRIVATE multifloat-fortran)
add_executable(multifloat_test ./fsrc/multifloat_test.f90)
target_link_libraries(multifloat_test PRIVATE multifloat-fortran)
