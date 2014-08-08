#! /bin/bash

# Insert the following lines into your .bashrc
#   export ARMCI_HOME=/home/pablo/Applications/ga-v5.3/installTree
#   export PARMETIS_HOME=/usr/
#   export LD_LIBRARY_PATH=/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server:$ARMCI_HOME/lib:$LD_LIBRARY_PATH

export ROSE_HOME=/home/pablo/Applications/rose-v0.9.5a/installTree/
export BOOST_HOME=/home/pablo/Applications/boost-v1.47.0/installTree/
export MPI_HOME=/usr/bin/

cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=/usr $PWD
make
sudo make install
