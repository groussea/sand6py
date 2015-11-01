
include_directories( ../src )

add_executable ( d6 d6.cc )
target_link_libraries ( d6 d6_simu )

add_executable ( d62vtk d62vtk.cc )
target_link_libraries ( d62vtk d6_visu )

if( OpenGL_FOUND )
	## GL Viewer
	find_package( QGLVIEWER )
	if( QGLVIEWER_FOUND )
		find_package(Qt4 REQUIRED QtGui QtXml QtOpenGL )

		include_directories( ${QGLVIEWER_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
		include_directories( ${QT_QTGUI_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR}  ${QT_QTOPENGL_INCLUDE_DIR} )
		SET( GL_LIBRARIES ${QGLVIEWER_LIBRARIES} ${OPENGL_LIBRARIES} Qt4::QtGui Qt4::QtXml Qt4::QtOpenGL )

		add_executable ( d6gl d6gl.cc GLViewer.cc GLViewer.hh )
		target_link_libraries ( d6gl d6_visu d6_gl ${GL_LIBRARIES}  )
	else()
	message( STATUS "QGLViewer not found, not builder viewer" )
	endif()
endif()
