include(${CMAKE_SOURCE_DIR}/cmake/FindRunTimeDependencies.cmake)

set(CMAKE_CXX_COMPILER ${MPICXX})
set(CMAKE_C_COMPILER ${MPICC})

add_subdirectory(include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ARMCI_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIR})
include_directories(${PATOH_INCLUDE_DIR})

add_subdirectory(DataStructures)

set(RUNTIME_SRCS
	${CMAKE_CURRENT_SOURCE_DIR}/access_data.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cdefs.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/communicator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/external.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/global_data.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/global_loop.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hypergraph.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Inspector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inspector_metis.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/local_comm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/local_data.cpp
	)


add_library(${PROJECT_NAME} ${RUNTIME_SRCS} ${DATA_STRUCTURES_SRCS} ${RUNTIME_INCLUDES} ${DATA_STRUCTURES_INCLUDES})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION 0.1)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${PATOH_LIBS} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/ie.h DESTINATION include)
install(FILES ${PATOH_INCLUDE_DIR}/patoh.h DESTINATION include)
