
MESSAGE( STATUS "Updating git revision" )

SET( GIT_CMD    git )
SET( D6_VERSION "1.0-review" )

# Get the current working branch
execute_process(
  COMMAND ${GIT_CMD} rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  RESULT_VARIABLE GIT_OK
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if( NOT( ${GIT_OK} EQUAL 0 ) )
	SET( GIT_BRANCH "release" )
endif()

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND ${GIT_CMD} log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  RESULT_VARIABLE GIT_OK
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if( NOT( ${GIT_OK} EQUAL 0 ) )
	SET( GIT_COMMIT_HASH ${D6_VERSION} )
endif()

if( ${CMAKE_VERSION} VERSION_GREATER 3.0.0 )
  string(TIMESTAMP COMPILE_TIME "%Y-%m-%d %H:%M")
endif()

configure_file(${TEMPLATE} ${OUTPUT})

