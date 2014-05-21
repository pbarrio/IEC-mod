include(${CMAKE_SOURCE_DIR}/cmake/FindCompileTimeDependencies.cmake)
find_program(GPP44 NAMES g++-4.4 g++ PATHS /usr/bin /usr/local/bin $ENV{GPP44_DIR} )
message("GPP44:"${GPP44})
set(CMAKE_CXX_COMPILER ${GPP44})

add_subdirectory(include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(ScopExtractor)

set(COMPILETIME_SRCS
	${CMAKE_CURRENT_SOURCE_DIR}/access_details.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/analyseAST.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/array_details.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/aux.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/codegen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/conditional.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/data_scalars.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/detector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/driver.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/executor_gen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/indarray_details.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inspector_gen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/loops.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/safe_scalars.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/unit_stride.cpp
	)

add_executable(${PROJECT_NAME}c ${SCOPEXTRACTOR_SRCS} ${COMPILETIME_SRCS} ${SCOPEXTRACTOR_INCLUDES} ${COMPILETIME_INCLUDES} )
target_link_libraries(${PROJECT_NAME}c ${ROSE_LIBS} ${BOOST_LIBS})
set_target_properties(${PROJECT_NAME}c PROPERTIES VERSION 0.1)
install(TARGETS ${PROJECT_NAME}c DESTINATION bin)
