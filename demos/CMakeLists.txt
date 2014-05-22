cmake_minimum_required(VERSION 2.8)
include(${CMAKE_SOURCE_DIR}/../cmake/FindRunTimeDependencies.cmake)

#IE LIBRARY and INCLUDES
find_library(IE_LIBS ie PATH ${IE_HOME} $ENV{IE_HOME} ${IE_HOME}/lib $ENV{IE_HOME}/lib)
find_path(IE_INCLUDE_DIR ie.h PATH ${IE_HOME} $ENV{IE_HOME} ${IE_HOME}/include $ENV{IE_HOME}/include)
find_program(IE_COMPILER NAME iec PATHS ${IE_HOME} ${IE_HOME}/bin $ENV{IE_HOME} $ENV{IE_HOME}/bin)

include_directories(${IE_INCLUDE_DIR})

message("IE : "${IE_LIBS},${IE_INCLUDE_DIR})

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
set(CMAKE_C_LINK_EXECUTABLE "<CMAKE_CXX_COMPILER> <FLAGS> <CMAKE_C_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")

add_subdirectory(smvp)
add_subdirectory(quake)
add_subdirectory(cg)
add_subdirectory(p3)
