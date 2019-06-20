#!/bin/bash

SEQ=list
OUT=time.csv

rm $OUT

while read NP ; do

	LINE=$(cat out/$NP | grep '^stats ' | sed 's/[^ ]*=//g;s/stats //g' | tail -1)

	echo $(($NP*$NP)) $LINE >> $OUT

done < $SEQ
