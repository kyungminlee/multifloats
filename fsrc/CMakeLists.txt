project(multifloatsf LANGUAGES Fortran)

include(GNUInstallDirs)
include(FortranCompiler)

# Detect usable REAL(16) support. The test exercises real(v, 16),
# ieee_is_finite, and arithmetic — flang may accept the type declaration
# but reject intrinsics operating on it.
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
set(_FYPP_DEFS "")
if(NOT FORTRAN_HAS_REAL16)
    list(APPEND _FYPP_DEFS "-DMULTIFLOATS_NO_QUAD")
endif()

set(_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(_GENERATED_SRC ${_GENERATED_DIR}/multifloats.f90)

add_custom_command(
    OUTPUT ${_GENERATED_SRC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE} ${_FYPP_DEFS}
        ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
        ${_GENERATED_SRC}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/multifloats.fypp
    VERBATIM)

add_library(multifloatsf STATIC ${_GENERATED_SRC})
fortran_module_layout(multifloatsf)
target_include_directories(multifloatsf PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# Fat LTO for GCC; plain -flto for LLVM.
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multifloatsf PRIVATE -O3 -flto -ffat-lto-objects)
        target_link_options(multifloatsf PUBLIC -flto)
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
        target_compile_options(multifloatsf PRIVATE -O3 -flto)
        target_link_options(multifloatsf PUBLIC -flto)
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        # ifx defaults to -fp-model=fast at -O2+. DD error-free
        # transformations require strict IEEE associativity and no FMA
        # contraction — see src/CMakeLists.txt for the full rationale.
        target_compile_options(multifloatsf PRIVATE -O3 -fp-model=precise)
    endif()
else()
    target_compile_options(multifloatsf PRIVATE -O2)
endif()

# Install: compiler-tagged library + .mod files + Config.cmake
fortran_install_library(multifloatsf)
fortran_install_modules(multifloatsf)

# Also install the generated Fortran source alongside the installed .mod
# files so downstream projects using a mismatched compiler can recompile
# the module locally. .f90 source is compiler-independent, so it is
# placed at <libdir>/fmod/ (no compiler-version subdirectory).
install(FILES ${_GENERATED_SRC}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/fmod)
