#! /bin/bash

export ROSE_HOME=/home/pbarrio/apps/rose-v0.9.5a/installTree/
export BOOST_HOME=/usr/
export MPI_HOME=/usr/local/bin/
export ARMCI_HOME=/home/pbarrio/apps/ga-5-3/installDir/
export PARMETIS_HOME=/home/pbarrio/apps/parmetis-4.0.3/installDir/

export LD_LIBRARY_PATH=/usr/lib/jvm/java-6-openjdk-amd64/jre/lib/amd64/server:$LD_LIBRARY_PATH

cmake -DCMAKE_INSTALL_PREFIX=/usr/local $PWD
make install
