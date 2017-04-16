#! /bin/bash

NPROC=15

if [ $# != 1 ]; then
	echo "ERROR: must specify a process ID for perf"
	exit -1
fi

np1=$(($1))
np2=$(($NPROC - 1 - $1))

if [ $np1 -lt 0 ] || [ $np2 -lt 0 ]; then
	echo "ERROR: the process must be within range [0, $(($NPROC - 1))]"
	exit -1
fi

if [ $np1 -gt 0 ]; then
	cmd1="-np $np1 iec_quake_block pabloutp.in :"
fi
if [ $np2 -gt 0 ]; then
	cmd2=": -np $np2 iec_quake_block pabloutp.in"
fi

echo "Gathering performance of process $1"
mkdir -p results
mpirun \
	$cmd1 \
	-np 1 perf record -o results/p$1.perf ./iec_quake_block pabloutp.in \
	$cmd2
