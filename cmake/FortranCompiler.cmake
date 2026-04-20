# FortranCompiler.cmake
#
# Detects the Fortran compiler's .mod file format version and provides
# functions for compiler-aware installation of Fortran modules and libraries.
#
# After include(FortranCompiler), the following variables are set:
#
#   FORTRAN_COMPILER_FAMILY   - Normalized compiler family (gfortran, intel, flang, nvhpc, nag, cray)
#   FORTRAN_COMPILER_VERSION  - Full compiler version string
#   FORTRAN_MOD_VERSION       - Internal .mod format version integer (or "unknown")
#   FORTRAN_MOD_COMPAT_TAG    - Module compat tag for .mod dirs (e.g. gfortran-mod15)
#   FORTRAN_COMPILER_TAG      - Compiler version tag for libraries (e.g. gfortran-14)
#
# Functions provided:
#
#   fortran_module_layout(<target>)
#     Sets up build-time and install-time module directories for the target.
#     Must be called before fortran_install_modules() or fortran_install_library().
#
#   fortran_install_modules(<target> [DESTINATION <base>])
#     Installs .mod/.smod files to <base>/fmod/<mod-tag>/.
#     Requires fortran_module_layout() to have been called on the target first.
#
#   fortran_install_library(<target> [NAMESPACE <ns>] [EXPORT <export-name>])
#     Installs the library with a compiler-tagged filename and generates a
#     Config.cmake for find_package() support.

if(_FORTRAN_COMPILER_INCLUDED)
  return()
endif()
set(_FORTRAN_COMPILER_INCLUDED TRUE)

# ---------------------------------------------------------------------------
# Verify Fortran is enabled
# ---------------------------------------------------------------------------
get_property(_fc_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if(NOT "Fortran" IN_LIST _fc_languages)
  message(FATAL_ERROR "FortranCompiler: Fortran language must be enabled before including this module.")
endif()
unset(_fc_languages)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Detect compiler family and ABI version tag
#
# The snippet below reads CMAKE_Fortran_COMPILER_ID / *_VERSION and sets
# _fc_family and _fc_abi_version. It is used here at configure time (for
# this build's producer tag) and embedded verbatim into the generated
# Config.cmake below (for the consumer's tag at find_package() time). The
# shared source keeps producer- and consumer-side family/version mapping
# from drifting.
#
# Compiler IDs used by CMake:
#   GNU       - gfortran
#   Intel     - ifort (classic)
#   IntelLLVM - ifx
#   LLVMFlang - LLVM Flang (flang-new, the official LLVM Fortran compiler)
#   Flang     - Classic Flang (PGI-derived, incompatible with LLVM Flang)
#   NVHPC     - NVIDIA nvfortran (PGI lineage, shares classic Flang .mod format)
#   NAG       - NAG Fortran
#   Cray      - Cray Fortran (CCE)
# ---------------------------------------------------------------------------
set(_FC_DETECT_CODE [[
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(_fc_family "gfortran")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(_fc_family "intel")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(_fc_family "intel")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(_fc_family "flang")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
  set(_fc_family "flang-classic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(_fc_family "nvhpc")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set(_fc_family "nag")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(_fc_family "cray")
else()
  set(_fc_family "${CMAKE_Fortran_COMPILER_ID}")
  string(TOLOWER "${_fc_family}" _fc_family)
endif()

# Version truncation per family:
#   gfortran / flang -> major only (ABI stable within a release series)
#   intel            -> major.minor (ABI can change at minor releases)
#   others           -> full version (conservative)
if(_fc_family STREQUAL "gfortran" OR _fc_family STREQUAL "flang")
  string(REGEX MATCH "^([0-9]+)" _fc_abi_version "${CMAKE_Fortran_COMPILER_VERSION}")
elseif(_fc_family STREQUAL "intel")
  string(REGEX MATCH "^([0-9]+\\.[0-9]+)" _fc_abi_version "${CMAKE_Fortran_COMPILER_VERSION}")
else()
  set(_fc_abi_version "${CMAKE_Fortran_COMPILER_VERSION}")
endif()
]])

set(FORTRAN_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
cmake_language(EVAL CODE "${_FC_DETECT_CODE}")
set(FORTRAN_COMPILER_FAMILY "${_fc_family}")
unset(_fc_family)

# ---------------------------------------------------------------------------
# Determine .mod format version from compiler family + version
# ---------------------------------------------------------------------------
set(FORTRAN_MOD_VERSION "unknown")

if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran")
  # GCC major version -> MOD_VERSION
  # GCC < 4.4: unversioned .mod format (unsupported)
  string(REGEX MATCH "^([0-9]+)" _fc_gcc_major "${FORTRAN_COMPILER_VERSION}")
  if(_fc_gcc_major VERSION_GREATER_EQUAL 15)
    set(FORTRAN_MOD_VERSION "16")
  elseif(_fc_gcc_major VERSION_GREATER_EQUAL 8)
    set(FORTRAN_MOD_VERSION "15")
  elseif(_fc_gcc_major VERSION_GREATER_EQUAL 5)
    set(FORTRAN_MOD_VERSION "14")
  elseif(_fc_gcc_major EQUAL 4)
    string(REGEX MATCH "^4\\.([0-9]+)" _fc_gcc_4minor "${FORTRAN_COMPILER_VERSION}")
    set(_fc_minor "${CMAKE_MATCH_1}")
    if(_fc_minor EQUAL 9)
      set(FORTRAN_MOD_VERSION "12")
    elseif(_fc_minor EQUAL 8)
      set(FORTRAN_MOD_VERSION "10")
    elseif(_fc_minor EQUAL 7)
      set(FORTRAN_MOD_VERSION "9")
    elseif(_fc_minor EQUAL 6)
      set(FORTRAN_MOD_VERSION "6")
    elseif(_fc_minor EQUAL 5)
      set(FORTRAN_MOD_VERSION "4")
    elseif(_fc_minor EQUAL 4)
      set(FORTRAN_MOD_VERSION "0")
    endif()
    unset(_fc_minor)
    unset(_fc_gcc_4minor)
  endif()
  unset(_fc_gcc_major)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Classic ifort versioning
  # Version numbering jumped from 19.x to 2021.x (year-based oneAPI scheme).
  # All versions 18.x through 2021.x compare correctly with VERSION_GREATER_EQUAL
  # because 18 < 2020 < 2021 numerically.
  if(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2021.10")
    set(FORTRAN_MOD_VERSION "13")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "18.0")
    set(FORTRAN_MOD_VERSION "12")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "17.0")
    set(FORTRAN_MOD_VERSION "11")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "16.0")
    set(FORTRAN_MOD_VERSION "10")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # ifx versioning (also uses year-based oneAPI scheme)
  if(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2023.2")
    set(FORTRAN_MOD_VERSION "13")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2021.0")
    set(FORTRAN_MOD_VERSION "12")
  endif()

elseif(FORTRAN_COMPILER_FAMILY STREQUAL "flang")
  # LLVM Flang: text-based .mod files with !mod$ v1 header
  set(FORTRAN_MOD_VERSION "1")

# flang-classic, nvhpc, nag, cray: no known .mod version mapping.
# They fall through to FORTRAN_MOD_VERSION = "unknown" and get tagged
# by full compiler version (conservative but safe).
endif()

# ---------------------------------------------------------------------------
# Build tags:
#   FORTRAN_MOD_COMPAT_TAG  - for .mod directories (by format compatibility)
#   FORTRAN_COMPILER_TAG    - for library files (by ABI-relevant version)
# ---------------------------------------------------------------------------
set(FORTRAN_COMPILER_TAG "${FORTRAN_COMPILER_FAMILY}-${_fc_abi_version}")
unset(_fc_abi_version)

if(FORTRAN_MOD_VERSION STREQUAL "unknown")
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_TAG}")
else()
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_FAMILY}-mod${FORTRAN_MOD_VERSION}")
endif()

message(STATUS "FortranCompiler: compiler=${CMAKE_Fortran_COMPILER_ID} ${FORTRAN_COMPILER_VERSION}")
message(STATUS "FortranCompiler: family=${FORTRAN_COMPILER_FAMILY}, mod_version=${FORTRAN_MOD_VERSION}")
message(STATUS "FortranCompiler: mod_tag=${FORTRAN_MOD_COMPAT_TAG}, lib_tag=${FORTRAN_COMPILER_TAG}")

# ---------------------------------------------------------------------------
# Helper: ensure INSTALL_INTERFACE paths are relative (required for
# relocatable packages). GNUInstallDirs can produce absolute paths on
# some platforms.
# ---------------------------------------------------------------------------
function(_fc_make_relative_installdir outvar path)
  if(IS_ABSOLUTE "${path}")
    file(RELATIVE_PATH _fc_relpath "${CMAKE_INSTALL_PREFIX}" "${path}")
    set(${outvar} "${_fc_relpath}" PARENT_SCOPE)
  else()
    set(${outvar} "${path}" PARENT_SCOPE)
  endif()
endfunction()

# ---------------------------------------------------------------------------
# fortran_module_layout(<target>)
#
# Configures the target's module output directory and include paths.
# Uses FORTRAN_MOD_COMPAT_TAG for module directories.
# Must be called before fortran_install_modules() or fortran_install_library().
# ---------------------------------------------------------------------------
function(fortran_module_layout target)
  set(_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}")

  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY "${_moddir}"
  )

  # Ensure the install interface path is relative for relocatable packages
  _fc_make_relative_installdir(_fc_rel_libdir "${CMAKE_INSTALL_LIBDIR}")

  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${_moddir}>
    $<INSTALL_INTERFACE:${_fc_rel_libdir}/fmod/${FORTRAN_MOD_COMPAT_TAG}>
  )
endfunction()

# ---------------------------------------------------------------------------
# fortran_install_modules(<target> [DESTINATION <base>])
#
# Installs .mod and .smod files to <base>/fmod/<mod-tag>/
# Default DESTINATION is ${CMAKE_INSTALL_LIBDIR}.
# Requires fortran_module_layout() to have been called on the target first.
# ---------------------------------------------------------------------------
function(fortran_install_modules target)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "DESTINATION" "")
  if(ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "fortran_install_modules: unexpected arguments: ${ARG_UNPARSED_ARGUMENTS}")
  endif()

  # Validate that fortran_module_layout was called
  get_target_property(_fc_moddir ${target} Fortran_MODULE_DIRECTORY)
  if(NOT _fc_moddir)
    message(FATAL_ERROR
      "fortran_install_modules(${target}): Fortran_MODULE_DIRECTORY is not set. "
      "Call fortran_module_layout(${target}) first.")
  endif()

  if(NOT ARG_DESTINATION)
    set(ARG_DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  endif()

  set(_install_moddir "${ARG_DESTINATION}/fmod/${FORTRAN_MOD_COMPAT_TAG}")

  install(
    DIRECTORY "${_fc_moddir}/"
    DESTINATION "${_install_moddir}"
    FILES_MATCHING
      PATTERN "*.mod"
      PATTERN "*.smod"
  )
endfunction()

# ---------------------------------------------------------------------------
# fortran_install_library(<target>
#     [NAMESPACE <ns>]
#     [EXPORT <export-name>]
#     [DESTINATION <lib-dir>])
#
# Installs the library with a compiler-version-tagged filename (for ABI),
# while the export/config system uses the mod compat tag to find modules.
#
# Note: This function generates a <ProjectName>Config.cmake. If your project
# has multiple Fortran library targets, add them all to a single EXPORT set
# by passing the same EXPORT name. The Config.cmake and export file are
# generated only once per export set.
# ---------------------------------------------------------------------------
function(fortran_install_library target)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "NAMESPACE;EXPORT;DESTINATION" "")
  if(ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "fortran_install_library: unexpected arguments: ${ARG_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARG_NAMESPACE)
    set(ARG_NAMESPACE "${PROJECT_NAME}::")
  endif()
  if(NOT ARG_EXPORT)
    set(ARG_EXPORT "${PROJECT_NAME}Targets")
  endif()
  if(NOT ARG_DESTINATION)
    set(ARG_DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  endif()

  set(_config_name "${PROJECT_NAME}")
  set(_cmake_install_dir "${ARG_DESTINATION}/cmake/${_config_name}")
  set(_targets_file "${ARG_EXPORT}-${FORTRAN_COMPILER_TAG}.cmake")

  # Tag the library output filename by compiler version (ABI compatibility)
  set_target_properties(${target} PROPERTIES
    OUTPUT_NAME "${target}-${FORTRAN_COMPILER_TAG}"
  )

  # Add target to the export set
  install(TARGETS ${target}
    EXPORT ${ARG_EXPORT}
    ARCHIVE DESTINATION "${ARG_DESTINATION}"
    LIBRARY DESTINATION "${ARG_DESTINATION}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )

  # Only generate the export file and Config.cmake once per export set.
  # Multiple targets can share the same EXPORT name; CMake accumulates them
  # from all install(TARGETS ... EXPORT <name>) calls. But install(EXPORT)
  # and Config.cmake generation must happen exactly once.
  get_property(_fc_config_generated GLOBAL PROPERTY _FC_CONFIG_GENERATED_${ARG_EXPORT})
  if(_fc_config_generated)
    return()
  endif()
  set_property(GLOBAL PROPERTY _FC_CONFIG_GENERATED_${ARG_EXPORT} TRUE)

  # Write the export set (includes all targets added to this EXPORT name)
  install(EXPORT ${ARG_EXPORT}
    FILE "${_targets_file}"
    NAMESPACE "${ARG_NAMESPACE}"
    DESTINATION "${_cmake_install_dir}"
  )

  # Generate Config.cmake that finds the right targets file.
  # Strategy: derive the consumer's compiler tag and look for an exact match.
  # If no exact match, fail with an informative error listing available builds.
  set(_config_content "\
# ${_config_name}Config.cmake
# Auto-generated by FortranCompiler.cmake
#
# Detects the consuming compiler and includes the matching targets file.
# Library files are tagged by compiler version (ABI compatibility).
# Module directories are tagged by .mod format version (compile-time compatibility).

cmake_minimum_required(VERSION 3.12)

# --- Derive consumer's compiler family and ABI version tag ---
# Detection snippet shared with FortranCompiler.cmake at producer build
# time — see _FC_DETECT_CODE there. Sets _fc_family, _fc_abi_version from
# CMAKE_Fortran_COMPILER_ID / *_VERSION of the consumer.
${_FC_DETECT_CODE}
set(_FC_consumer_tag \"\${_fc_family}-\${_fc_abi_version}\")
unset(_fc_family)
unset(_fc_abi_version)

# Look for exact compiler version match
set(_FC_targets_file \"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}-\${_FC_consumer_tag}.cmake\")

if(NOT EXISTS \"\${_FC_targets_file}\")
  # No exact match — list available builds and fail
  file(GLOB _FC_available \"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}-*.cmake\")
  set(_FC_available_names \"\")
  foreach(_FC_f IN LISTS _FC_available)
    get_filename_component(_FC_fname \"\${_FC_f}\" NAME)
    list(APPEND _FC_available_names \"\${_FC_fname}\")
  endforeach()
  list(JOIN _FC_available_names \", \" _FC_available_list)
  set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    \"${_config_name}: no pre-built library found for compiler '\${_FC_consumer_tag}'. Available: [\${_FC_available_list}]\")
  unset(_FC_consumer_family)
  unset(_FC_consumer_tag)
  unset(_FC_targets_file)
  unset(_FC_available)
  unset(_FC_available_names)
  unset(_FC_available_list)
  return()
endif()

include(\"\${_FC_targets_file}\")

unset(_FC_consumer_family)
unset(_FC_consumer_tag)
unset(_FC_targets_file)
")

  file(GENERATE
    OUTPUT "${PROJECT_BINARY_DIR}/cmake/${_config_name}Config.cmake"
    CONTENT "${_config_content}"
  )

  install(
    FILES "${PROJECT_BINARY_DIR}/cmake/${_config_name}Config.cmake"
    DESTINATION "${_cmake_install_dir}"
  )
endfunction()
