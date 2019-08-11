MACRO(ADD_D6_LIB lib )

  execute_process( COMMAND
	${CMAKE_COMMAND}
	-DD6_LIB=${lib}
	-DD6_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
	-DD6_DIM_STR=${DIM_STR}
	-P ${CMAKE_SOURCE_DIR}/cmake/UpdateSources.cmake )

  file(STRINGS ${lib}/CMakeSources_${DIM_STR}.txt SRC_LIST)

  add_library(${CMAKE_PROJECT_NAME}_${lib} SHARED ${SRC_LIST})
  target_include_directories(${CMAKE_PROJECT_NAME}_${lib} BEFORE
	  PUBLIC ${DIM_STR}/${lib} ${lib} )
  target_link_libraries(${CMAKE_PROJECT_NAME}_${lib} ${ARGN} )

ENDMACRO( ADD_D6_LIB )

INCLUDE_DIRECTORIES( . ${DIM_STR} )

ADD_D6_LIB(utils)
ADD_D6_LIB(geo  d6_utils ${EXTERNAL_LIBRARIES}  )
ADD_D6_LIB(solve ${EXTERNAL_LIBRARIES} )
ADD_D6_LIB(simu d6_utils d6_geo d6_solve ${EXTERNAL_LIBRARIES} )

ADD_D6_LIB(mono d6_simu ${EXTERNAL_LIBRARIES} )
ADD_D6_LIB(diphasic d6_mono ${EXTERNAL_LIBRARIES} )

ADD_D6_LIB(visu d6_utils d6_geo ${EXTERNAL_LIBRARIES} )
if( OPENGL_FOUND )
	ADD_D6_LIB(gl d6_utils d6_geo d6_visu ${OPENGL_LIBRARIES})
endif()
