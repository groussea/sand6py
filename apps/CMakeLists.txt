
include_directories( ../src )

add_executable ( d6 d6.cc )
target_link_libraries ( d6 d6_utils d6_simu )

add_executable ( d62vtk d62vtk.cc )
target_link_libraries ( d62vtk d6_utils d6_visu )


