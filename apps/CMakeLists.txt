

# Git info
SET( D6_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/generated/d6_version.cc )

add_custom_target(d6_version ALL
	${CMAKE_COMMAND}
	-DTEMPLATE=${CMAKE_CURRENT_SOURCE_DIR}/d6_version.tpl.cc
	-DOUTPUT=${D6_VERSION_FILE}
	-P ${CMAKE_SOURCE_DIR}/cmake/GenVersion.cmake )

set_source_files_properties(${D6_VERSION_FILE}
	PROPERTIES GENERATED TRUE )

# Executables

include_directories( ../src ../src/${DIM_STR} )

add_executable ( d6 d6.cc ${D6_VERSION_FILE} )
add_dependencies( d6 d6_version )
target_link_libraries ( d6 d6_mono )

add_executable ( dd6 dd6.cc ${D6_VERSION_FILE} )
add_dependencies( dd6 d6_version )
target_link_libraries ( dd6 d6_diphasic )

add_executable ( d62vtk d62vtk.cc )
target_link_libraries ( d62vtk d6_visu )

add_executable ( d6_solvePrimal d6_solvePrimal.cc )
target_link_libraries ( d6_solvePrimal d6_utils d6_solve )

add_executable ( dd6_solvePrimal dd6_solvePrimal.cc )
target_link_libraries ( dd6_solvePrimal d6_utils d6_diphasic )

if(DIM2)
	add_executable ( dd6_plot 2d/dd6_plot.cc )
	target_link_libraries ( dd6_plot d6_visu )

	if( OPENGL_FOUND )
		## GL Viewer

		## GLFW
		find_package(PkgConfig)
		if( PkgConfig_FOUND )
			pkg_search_module(GLFW glfw3)
		endif()
		if( GLFW_FOUND )
		  add_definitions(-DGLFW3)
	    else()
		  find_package(GLFW)
	    endif()

		if( GLFW_FOUND )
			MESSAGE( STATUS "Found glfw" )
			include_directories(${GLFW_INCLUDE_DIRS})
			link_directories(${GLFW_LIBRARY_DIRS})
			message(STATUS ${GLFW_LIBRARIES})
			SET( D6_GL_LIBRARIES ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES} )
		else()
			message( STATUS "GLFW not found, not building viewer" )
		endif()

	endif()

else( DIM2 )
	add_executable ( d6_analyze 3d/d6_analyze.cc )
	target_link_libraries ( d6_analyze d6_visu )

	#if( OPENGL_FOUND )
	#	## GL Viewer
	#	find_package( QGLVIEWER )
	#	if( QGLVIEWER_FOUND )
	#		find_package(Qt4 REQUIRED QtGui QtXml QtOpenGL )

	#		include_directories( ${QGLVIEWER_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
	#		include_directories( ${QT_QTGUI_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR}  ${QT_QTOPENGL_INCLUDE_DIR} )
	#		SET( D6_GL_LIBRARIES ${QGLVIEWER_LIBRARIES} ${OPENGL_LIBRARIES} Qt4::QtGui Qt4::QtXml Qt4::QtOpenGL )

	#	else()
	#		message( STATUS "QGLViewer not found, not building viewer" )
	#	endif()
	#endif()
		## GLFW
		find_package(PkgConfig)
		if( PkgConfig_FOUND )
			pkg_search_module(GLFW glfw3)
		endif()
		if( GLFW_FOUND )
		  add_definitions(-DGLFW3)
	    else()
		  find_package(GLFW)
	    endif()

		if( GLFW_FOUND )
			MESSAGE( STATUS "Found glfw" )
			include_directories(${GLFW_INCLUDE_DIRS})
			link_directories(${GLFW_LIBRARY_DIRS})
			message(STATUS ${GLFW_LIBRARIES})
			SET( D6_GL_LIBRARIES ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES} )
		else()
			message( STATUS "GLFW not found, not building viewer" )
		endif()
endif( DIM2 )

if( D6_GL_LIBRARIES )

if(DIM2)
	add_executable ( d6gl_2d ${DIM_STR}/d6gl.cc ${DIM_STR}/GLViewer.cc ${DIM_STR}/GLViewer.hh )
	target_link_libraries ( d6gl_2d d6_visu d6_gl ${D6_GL_LIBRARIES}  )
else( DIM2 )
	add_executable ( d6gl ${DIM_STR}/d6gl.cc ${DIM_STR}/GLViewer.cc ${DIM_STR}/GLViewer.hh )
		target_link_libraries ( d6gl d6_visu d6_gl ${D6_GL_LIBRARIES}  )
endif( DIM2 )	
endif()

