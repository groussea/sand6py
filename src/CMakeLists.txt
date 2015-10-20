MACRO(ADD_D6_LIB lib )

  execute_process( COMMAND ../updateCMakeSources.rb ${lib} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  file(STRINGS ${lib}/CMakeSources.txt SRC_LIST)

  include_directories(. ${lib})
  add_library(${CMAKE_PROJECT_NAME}_${lib} SHARED ${SRC_LIST})

ENDMACRO( ADD_D6_LIB )

ADD_D6_LIB(utils)
ADD_D6_LIB(geo)
#ADD_D6_LIB(solve)
ADD_D6_LIB(simu)
#ADD_D6_LIB(visu)

