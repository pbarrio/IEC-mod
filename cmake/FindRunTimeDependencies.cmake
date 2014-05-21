#Find MPI
find_program(MPICXX mpicxx)
find_program(MPICC mpicc)
find_program(MPIF90 mpif90)
find_path(MPI_INCLUDE_DIR mpi.h PATH ${MPI_HOME}/include $ENV{MPI_HOME}/include)
message("Compiler : "${MPICXX},${MPICC},${MPIF90})

#Find ARMCI/GA
find_library(ARMCI_LIBS armci PATH ${ARMCI_HOME}/lib $ENV{ARMCI_HOME}/lib)
find_path(ARMCI_INCLUDE_DIR armci.h ${ARMCI_HOME}/include $ENV{ARMCI_HOME}/include)
message("ARMCI Lib/Include : "${ARMCI_LIBS},${ARMCI_INCLUDE_DIR})

#Find Parmetis
find_library(PARMETIS_LIBS parmetis PATH ${PARMETIS_HOME} ${PARMETIS_HOME}/lib $ENV{PARMETIS_HOME} $ENV{PARMETIS_HOME}/lib)
find_library(METIS_LIBS metis PATH ${PARMETIS_HOME} ${PARMETIS_HOME}/lib $ENV{PARMETIS_HOME} $ENV{PARMETIS_HOME}/lib)
find_path(PARMETIS_INCLUDE_DIR parmetis.h PATH ${PARMETIS_HOME} ${PARMETIS_HOME}/include $ENV{PARMETIS_HOME} $ENV{PARMETIS_HOME}/include)
message("PARMETIS Lib/Include : "${PARMETIS_LIBS},${METIS_LIBS},${PARMETIS_INCLUDE_DIR})

#Find Patoh
find_library(PATOH_LIBS patoh PATH ${CMAKE_SOURCE_DIR}/external/PaToH ${IE_HOME}/lib $ENV{IE_HOME}/lib)
find_path(PATOH_INCLUDE_DIR patoh.h ${CMAKE_SOURCE_DIR}/external/PaToH ${IE_HOME}/include $ENV{IE_HOME}/include)

message("PATOH Lib/Include : "${PATOH_LIBS},${PATOH_INCLUDE_DIR})

#Find OpenMP
find_package(OpenMP)
message("OPENMP_FOUND ":${OPENMP_FOUND})
if( OPENMP_FOUND )
    message("OpenMP_C_FLAGS : "${OpenMP_C_FLAGS})
    message("OpenMP_CXX_FLAGS : "${OpenMP_CXX_FLAGS})
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif( OPENMP_FOUND )
