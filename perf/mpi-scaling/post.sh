#!/bin/bash

SEQ=list
OUT=time.csv

rm $OUT

while read NCPUS ; do

	LINE=$(cat out/$NCPUS | grep '^stats ' | sed 's/[^ ]*=//g;s/stats //g' | tail -1)

	NP=$((32/$NCPUS))
	RATIO=$(echo $NCPUS/32 | bc -l)
	PART=$(echo $LINE | awk '{printf "%e\n", $3 - $8}')

	echo $NCPUS $NP $RATIO $LINE $PART >> $OUT

done < $SEQ
