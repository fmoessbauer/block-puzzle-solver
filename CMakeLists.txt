cmake_minimum_required(VERSION 3.2)

project(blocks-puzzle)
set(SOURCES "main")

if(NOT "${TASK}" STREQUAL "")
  set(TASK_FLAG "-D${TASK} ")
endif()

set(RELEASE_COMPILE_FLAGS "-O3 -Ofast -Wall")
set(DEBUG_COMPILE_FLAGS "-g -Wall -Werr") 

add_executable(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${TASK_FLAG} ${DEBUG_COMPILE_FLAGS}")
else()
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${TASK_FLAG} ${RELEASE_COMPILE_FLAGS}")
endif()

if("$ENV{BLAZE_DIR}" STREQUAL "")
  set(BLAZE_DIR "$ENV{HOME}/opt/blaze-3.2")
else()
  set(BLAZE_DIR "$ENV{BLAZE_DIR}")
endif()

find_package( blaze REQUIRED) 
if( blaze_FOUND )
  add_library(blaze INTERFACE)
  target_link_libraries(blaze INTERFACE blaze::blaze)
endif()

find_package( Boost REQUIRED COMPONENTS graph)
if( Boost_FOUND )
  add_library(Boost INTERFACE)
  target_link_libraries(Boost INTERFACE Boost::graph)
endif()

target_link_libraries(${PROJECT_NAME} blaze Boost)

