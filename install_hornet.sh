#! /bin/bash

# Insert the following lines into your .bashrc
#   export ROSE_HOME=/home/pablo/Applications/rose-v0.9.5a/installTree/
export LD_LIBRARY_PATH=/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server:$LD_LIBRARY_PATH #:$ARMCI_HOME/lib:$ROSE_HOME/lib:$BOOST_HOME/lib

export MPI_HOME=/usr/

mkdir -p buildTree
mkdir -p installTree

cd buildTree
cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=../installTree .. \
    -DARMCI_HOME=/home/pablo/Applications/ga-v5.3/installTree \
    -DBOOST_HOME=/home/pablo/Applications/boost-v1.47.0/installTree/ \
    -DPARMETIS_HOME=/usr/ \
    -DROSE_HOME=/home/pablo/Applications/rose-v0.9.5a/installTree/
make -j4
make install
cd ..
