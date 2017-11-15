cmake_minimum_required(VERSION 3.2)

project(blocks-puzzle)
set(SOURCES "main")

add_executable(${PROJECT_NAME} ${SOURCES})
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)

if("$ENV{BLAZE_DIR}" STREQUAL "")
  set(BLAZE_DIR "$ENV{HOME}/opt/blaze-3.2")
else()
  set(BLAZE_DIR "$ENV{BLAZE_DIR}")
endif()

find_package( blaze REQUIRED HINTS ${BLAZE_DIR})                               

if( blaze_FOUND )
  message("blaze found")
  add_library(blaze INTERFACE)
  target_link_libraries(blaze INTERFACE blaze::blaze)
endif()
target_link_libraries(${PROJECT_NAME} blaze)

