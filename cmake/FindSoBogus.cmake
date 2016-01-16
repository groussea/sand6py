# - Try to find So-Bogus lib
#
#
# Once done this will define
#
#  SoBogus_FOUND - system has SoBogus lib
#  SoBogus_INCLUDE_DIR - the SoBogus include directory

find_path(SoBogus_INCLUDE_DIR NAMES bogus/Interfaces/Cadoux.hpp
  HINTS
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor
  PATHS
  $ENV{BOGUS_ROOT}
  ${CMAKE_INSTALL_PREFIX}/include
  )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SoBogus "SoBogus library not found ; set environment variable BOGUS_ROOT " SoBogus_INCLUDE_DIR)
mark_as_advanced( SoBogus_INCLUDE_DIR )
