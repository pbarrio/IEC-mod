#! /bin/bash

dir=`dirname $0`
echo $dir
for i in `seq 1 7000`; do
	echo "Start from $i"
	$dir/mesh-editor pabloinp.in delete nodes $i 7293 > test.in
	$dir/../orig_quake test.in > out 2>&1
	grep -n nan out
	if [ $? -eq 1 ]; then
		echo "Found a good mesh: size = $i - 1"
		break
	fi
done

rm test.in out
