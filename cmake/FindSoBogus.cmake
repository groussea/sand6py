# - Try to find So-Bogus lib
#
#
# Once done this will define
#
#  SoBogus_FOUND - system has SoBogus lib
#  SoBogus_INCLUDE_DIR - the SoBogus include directory
#  SoBogus_LIBRARIES - link this to use SoBogus

find_library(SoBogus_LIBRARIES NAMES bogus sobogus SoBogus bogus-static sobogus-static SoBogus-static
  HINTS 
  /usr/local/lib $ENV{BOGUS_ROOT}/lib
  PATH_SUFFIXES . Release Debug
  )
find_path(SoBogus_INCLUDE_DIR NAMES bogus/Interfaces/MecheInterface.hpp
  HINTS
  /usr/local/include $ENV{BOGUS_ROOT}/include
  )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SoBogus "SoBogus library not found ; set environment variable BOGUS_ROOT " SoBogus_LIBRARIES SoBogus_INCLUDE_DIR)
#make_library_set( ${SoBogus_LIBRARIES} )
mark_as_advanced( SoBogus_LIBRARIES )
mark_as_advanced( SoBogus_INCLUDE_DIR )

if( SOBOGUS_FOUND )
  if( ${SoBogus_LIBRARIES} MATCHES "-static\\." )
      SET( SOBOGUS_BUILT_AS_STATIC TRUE )
      add_definitions(-Dsobogus_BUILT_AS_STATIC )
  else()
      SET( SOBOGUS_BUILT_AS_STATIC FALSE )
  endif()

  string( REPLACE Release Debug SoBogus_LIB_DEBUG ${SoBogus_LIBRARIES} )
  string( REPLACE Debug Release SoBogus_LIB_RELEASE ${SoBogus_LIBRARIES} )

endif()
mark_as_advanced( SOBOGUS_BUILT_AS_STATIC )
