#! /bin/bash

# Insert the following lines into your .bashrc
#   export ROSE_HOME=/home/pablo/Applications/rose-v0.9.5a/installTree/
#   export LD_LIBRARY_PATH=/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server:$ARMCI_HOME/lib:$ROSE_HOME/lib:$BOOST_HOME/lib:$LD_LIBRARY_PATH

BUILD_DIR=$PWD/build
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BUILD_DIR=$PWD/build
INSTALL_DIR=$PWD/install

BOOST_HOME=/sw/openmpi/BOOST/1_49_0
MPI_HOME=/usr/
LD_LIBRARY_PATH=$BOOST_HOME/lib:$LD_LIBRARY_PATH

mkdir -p $BUILD_DIR
mkdir -p $INSTALL_DIR

cd $BUILD_DIR
cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR $ROOT_DIR \
    -DARMCI_HOME=/home/p097/PROJECT/lib/install-ga \
    -DBOOST_HOME=$BOOST_HOME \
    -DPARMETIS_HOME=/sw/openmpi/METIS/PARMETIS-4.0.3 \
    -DMETIS_HOME=/sw/openmpi/METIS/METIS-4.0.3 \
    -DROSE_HOME=/home/p097/PROJECT/lib/rose/installTree
make
make install
