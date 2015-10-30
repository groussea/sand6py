MACRO(ADD_D6_LIB lib )

  execute_process( COMMAND ../updateCMakeSources.rb ${lib} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  file(STRINGS ${lib}/CMakeSources.txt SRC_LIST)

  include_directories(. ${lib})
  add_library(${CMAKE_PROJECT_NAME}_${lib} SHARED ${SRC_LIST})
  target_link_libraries(${CMAKE_PROJECT_NAME}_${lib} ${ARGN} )

ENDMACRO( ADD_D6_LIB )


ADD_D6_LIB(utils)
ADD_D6_LIB(geo)
ADD_D6_LIB(solve ${EXTERNAL_LIBRARIES} )
ADD_D6_LIB(simu d6_utils d6_geo d6_solve ${EXTERNAL_LIBRARIES} )
ADD_D6_LIB(visu d6_utils d6_geo ${EXTERNAL_LIBRARIES} )

