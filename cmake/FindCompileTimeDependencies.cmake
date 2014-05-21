#Find Rose
find_library(ROSE_LIBS rose PATH ${ROSE_HOME}/lib $ENV{ROSE_HOME}/lib)
find_path(ROSE_INCLUDES rose.h PATH ${ROSE_HOME}/include $ENV{ROSE_HOME}/include)
include_directories(${ROSE_INCLUDES})

#Find Boost
find_library(BOOST_LIBS boost_regex PATH ${BOOST_HOME}/lib $ENV{BOOST_HOME}/lib)
find_path(BOOST_INCLUDES boost/regex.hpp PATH ${BOOST_HOME}/include $ENV{BOOST_HOME}/include)
include_directories(${BOOST_INCLUDES})

