
MESSAGE( STATUS "Updating sources for ${D6_LIB} ${D6_DIM_STR}" )

SET( GIT_CMD git )
SET( PATH ${D6_SRC_DIR}/${D6_LIB} )
SET( DEST ${PATH}/CMakeSources_${D6_DIM_STR}.txt )

execute_process(
  COMMAND ${GIT_CMD} ls-files ${PATH}  ${D6_SRC_DIR}/${D6_DIM_STR}/${D6_LIB}
  COMMAND grep -e "\\.\\(cc\\|hh\\|glsl\\)$"
  WORKING_DIRECTORY ${D6_SRC_DIR}
  RESULT_VARIABLE SOURCES_OK
  OUTPUT_VARIABLE SOURCES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(SOURCES_OK EQUAL 0)
    file( WRITE ${DEST} ${SOURCES})
endif()
