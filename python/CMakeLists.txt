
project(d6_python)
cmake_minimum_required(VERSION 2.8)
add_subdirectory(pybind11)
pybind11_add_module(d6_python src/d6_python.cpp)

include_directories( ../src ../src/${DIM_STR} )
#include_directories( ../vendor )
target_link_libraries(d6_python PRIVATE d6_mono)
target_link_libraries(d6_python PRIVATE d6_visu)
set_target_properties(d6_python PROPERTIES INSTALL_RPATH ${D6_BIN}/src)
set_target_properties(d6_python PROPERTIES SKIP_BUILD_RPATH FALSE)
# install(FILES d6py/__init__.py d6py/Tools.py DESTINATION . )
install(TARGETS d6_python LIBRARY DESTINATION . )

