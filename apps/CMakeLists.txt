

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

include_directories( ../src )

add_executable ( d6 d6.cc ${D6_VERSION_FILE} )
add_dependencies( d6 d6_version )
target_link_libraries ( d6 d6_simu )

add_executable ( d62vtk d62vtk.cc )
target_link_libraries ( d62vtk d6_visu )

if( OPENGL_FOUND )
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

