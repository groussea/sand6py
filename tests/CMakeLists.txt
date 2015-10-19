cmake_minimum_required(VERSION 2.8)

if( TESTS_NAMES )
    SET( SRC_LIST ${TESTS_NAMES} )
    MESSAGE( STATUS "Building tests: " ${SRC_LIST} )
else()
    execute_process( COMMAND ../updateCMakeSources.rb WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    file(STRINGS CMakeSources.txt SRC_LIST)
    MESSAGE( STATUS "Building all tests" )
endif()

include( ExternalProject )

if( MSVC AND MSVC_VERSION EQUAL 1700 )
  add_definitions( /D_VARIADIC_MAX=10 )
  set( GTEST_PATCH_COMMAND patch < ${CMAKE_CURRENT_SOURCE_DIR}/gtest-msvc11.patch )
  MESSAGE( STATUS "MSVC 2012 detected -- using gtest-msvc11.patch")
else()
  set( GTEST_PATCH_COMMAND "" )
endif()

if( MSVC )
    set( GTEST_CMAKE_ARGS -Dgtest_force_shared_crt=ON )
else()
    set( GTEST_CMAKE_ARGS )
endif()

# Add gtest
ExternalProject_Add(
    googletest
    URL https://github.com/google/googletest/archive/release-1.7.0.tar.gz
    URL_HASH MD5=4ff6353b2560df0afecfbda3b2763847
    CMAKE_ARGS ${GTEST_CMAKE_ARGS}
    UPDATE_COMMAND ""
    PATCH_COMMAND ${GTEST_PATCH_COMMAND}
    INSTALL_COMMAND ""
    )

SET( GTEST_LINK_LIBRARIES gtest gtest_main )
if(UNIX)
  SET( GTEST_LINK_LIBRARIES ${GTEST_LINK_LIBRARIES} pthread )
endif()

ExternalProject_Get_Property(googletest source_dir)
include_directories( ../src ${source_dir}/include)

ExternalProject_Get_Property(googletest binary_dir)
link_directories( ${binary_dir} )

add_executable(test${CMAKE_PROJECT_NAME} ${SRC_LIST})
add_test( ${CMAKE_PROJECT_NAME} test${CMAKE_PROJECT_NAME} )

add_dependencies(test${CMAKE_PROJECT_NAME} googletest)

target_link_libraries(test${CMAKE_PROJECT_NAME} ${EXTERNAL_LIBRARIES} ${GTEST_LINK_LIBRARIES} )
target_link_libraries(test${CMAKE_PROJECT_NAME} d6_utils d6_geo )

