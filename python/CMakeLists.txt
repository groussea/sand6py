

if (DIM2)
project(d6python2D)
else(DIM2)
project(d6python3D)
endif(DIM2)


# project(d6_python_2D)

cmake_minimum_required(VERSION 2.8)
add_subdirectory(pybind11)
if (DIM2)
pybind11_add_module(d6python_2D src/d6python2D.cpp)
include_directories( ../src ../src/${DIM_STR} )
#include_directories( ../vendor )
target_link_libraries(d6python2D PRIVATE d6_mono)
target_link_libraries(d6python2D PRIVATE d6_visu)
set_target_properties(d6python2D PROPERTIES INSTALL_RPATH ${D6_BIN}/src)
set_target_properties(d6python2D PROPERTIES SKIP_BUILD_RPATH FALSE)
# install(FILES d6py/__init__.py d6py/Tools.py DESTINATION . )
install(TARGETS d6python_2D LIBRARY DESTINATION . )
else(DIM2)
pybind11_add_module(d6python3D src/d6python3D.cpp)
include_directories( ../src ../src/${DIM_STR} )
#include_directories( ../vendor )
target_link_libraries(d6python3D PRIVATE d6_mono)
target_link_libraries(d6python3D PRIVATE d6_visu)
set_target_properties(d6python3D PROPERTIES INSTALL_RPATH ${D6_BIN}/src)
set_target_properties(d6python3D PROPERTIES SKIP_BUILD_RPATH FALSE)
# install(FILES d6py/__init__.py d6py/Tools.py DESTINATION . )
install(TARGETS d6python3D LIBRARY DESTINATION . )
endif(DIM2)




