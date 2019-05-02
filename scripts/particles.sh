#!/bin/bash

START=$1
STEP=$2
END=$3
CONF=$4

NP=$(seq $START $STEP $END)

for np in $NP; do
	echo "Running with ${np} particles"
	echo "@include \"$(basename $CONF)\"" > $CONF.$np
	echo "species = (
	{
	      particles = $np
	}
)" >> $CONF.$np
	./cpic -q $CONF.$np
done

