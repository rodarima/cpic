#!/bin/bash

START=$1
STEP=$2
END=$3
CONF=$4

BS=$(seq $START $STEP $END)

for bs in $BS; do
	echo "Running with ${bs}x${bs}"
	echo "@include \"$(basename $CONF)\"" > $CONF.$bs
	echo "grid = { blocksize = [$bs, $bs] }" >> $CONF.$bs
	./cpic -q $CONF.$bs
done

