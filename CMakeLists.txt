CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
set(PROJECT_NAME ie)
PROJECT(${PROJECT_NAME})
message("Project : "${PROJECT_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

add_subdirectory(src)
